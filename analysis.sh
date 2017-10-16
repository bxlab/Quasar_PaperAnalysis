#!/bin/bash
NCPUS=100
if [ "${NCPUS}" = "1" ]; then
  RUN="mpirun -np ${NCPUS}"
else
  RUN=""
fi

if [ ! -e Scripts/library.so ]; then
  cd Scripts && python setup.py build_ext --inplace && cd ../
fi

#rm Results/quasar_*_results.txt
mkdir -p Figures

# Generate Quality results
if [ ! -e Results/quasar_quality_results.txt ]; then
  echo "Quasar Quality"
  mkdir -p Quasar/Results
  A=($(eval grep -v dm6 Data/sample_key.txt | awk '{ print $1 }'))
  B=($(eval grep -v dm6 Data/sample_key.txt | awk '{ print $2 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${B[${I}]}
    ALIAS=${A[${I}]}
    C=(${ALIAS//_/ })
    if [ ! -e Data/${SAM}.hcd ]; then
      wget --no-check-certificate -O Data/${SAM}.hcd  https://bx.bio.jhu.edu/data/quasar/${SAM}.hcd
      wget --no-check-certificate -O Data/${SAM}_bin.hcp  https://bx.bio.jhu.edu/data/quasar/${SAM}_bin.hcp
      if [ ${SAM} != ${ALIAS} ]; then
        ln -s ${PWD}/Data/${SAM}_bin.hcp Data/${ALIAS}_bin.hcp
      fi
    fi
    if [ ! -e Data/${C[1]}_${C[3]}.fends ]; then
      wget --no-check-certificate -O Data/${C[1]}_${C[3]}.fends  https://bx.bio.jhu.edu/data/quasar/${C[1]}_${C[3]}.fends
    fi
    if [ ! -e Quasar/${ALIAS}.quasar ]; then
        echo "quality ${ALIAS}"
        ${RUN} hifive quasar -p Data/${ALIAS}_bin.hcp Quasar/${ALIAS}.quasar -r 1000000,200000,40000,10000 -d 0
    fi
    if [ ! -e Quasar/Results/${ALIAS}.txt ]; then
      echo "quality score ${ALIAS}"
      python Scripts/find_quasar_quality_score.py Quasar/${ALIAS}.quasar Quasar/Results/${ALIAS}.txt
    fi
  done

  A=($(eval grep dm6 Data/sample_key.txt | awk '{ print $1 }'))
  B=($(eval grep dm6 Data/sample_key.txt | awk '{ print $2 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${B[${I}]}
    ALIAS=${A[${I}]}
    C=(${SAM//_/ })
    if [ ! -e Data/${SAM}.hcd ]; then
      wget --no-check-certificate -O Data/${SAM}.hcd  https://bx.bio.jhu.edu/data/quasar/${SAM}.hcd
      wget --no-check-certificate -O Data/${SAM}_bin.hcp  https://bx.bio.jhu.edu/data/quasar/${SAM}_bin.hcp
      ln -s ${PWD}/Data/${SAM}_bin.hcp Data/${ALIAS}
    fi
    if [ ! -e Data/${C[1]}_${C[3]}.fends ]; then
      wget --no-check-certificate -O Data/${C[1]}_${C[3]}.fends  https://bx.bio.jhu.edu/data/quasar/${C[1]}_${C[3]}.fends
    fi
    if [ ! -e Quasar/${ALIAS}.quasar]; then
        echo "quality ${ALIAS}"
        ${RUN} hifive quasar -p Data/${ALIAS}_bin.hcp Quasar/${ALIAS}.quasar -r 100000,20000,4000 -d 0
    fi
    if [ ! -e Quasar/Results/${ALIAS}.txt ]; then
      echo "quality score ${ALIAS}"
      python Scripts/find_quasar_quality_score.py Quasar/${ALIAS}.quasar Quasar/Results/${ALIAS}.txt
    fi
  done

  A=($(eval cat Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  echo -e "Sample\tGenome\tCoverage\tResolution\tScore" > Results/quasar_quality_results.txt
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    B=(${SAM//_/ })
    GEN=${B[1]}
    if [ -e Quasar/Results/${SAM}.txt ]; then
      tail -n +2 Quasar/Results/${SAM}.txt | awk -v SAM="$SAM" -v GEN="$GEN" '{print SAM,GEN,$1,$2,$3}' OFS='\t' >> Results/quasar_quality_results.txt
    fi
  done
fi

# Generate Replicate results
if [ ! -e Results/quasar_replicate_results.txt ]; then
  echo "Quasar Replicate"
  A=($(eval grep -v dm6 Data/sample_key.txt | grep Rep1 | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]/_Rep1/}
    if [ ! -e Quasar/${SAM}_Pseudo.quasar ]; then
      echo "pseudo quality ${SAM}"
      ${RUN} python Scripts/QuasarPseudo.py Data/${SAM}_Rep1_bin.hcp Data/${SAM}_Rep2_bin.hcp Quasar/${SAM}_Pseudo.quasar -r 1000000,200000,40000,10000 -d 0
    fi
  done

  A=($(eval grep dm6 Data/sample_key.txt | grep Rep1 | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]/_Rep1/}
    if [ ! -e Quasar/${SAM}_Pseudo.quasar ]; then
      echo "pseudo quality ${SAM}"
      ${RUN} python Scripts/QuasarPseudo.py Data/${SAM}_Rep1_bin.hcp Data/${SAM}_Rep2_bin.hcp Quasar/${SAM}_Pseudo.quasar -r 100000,20000,4000 -d 0
    fi
  done

  echo -e "Sample1\tSample2\tGenome\tCoverage\tResolution\tScore" > Results/quasar_replicate_results.txt
  for I in hg38 mm10 dm6; do
    if [ ! -e Quasar/Results/quasar_replicate_results_${I}.txt ]; then
      echo "quality replicate ${I}"
      ${RUN} python Scripts/find_quasar_replicate_all_pair_scores.py Quasar/\*${I}\*.quasar Quasar/Results/quasar_replicate_results_${I}.txt
    fi
    tail -n +2 Quasar/Results/quasar_replicate_results_${I}.txt | awk -v GEN="$I" '{print $1,$2,GEN,$3,$4,$5}' OFS='\t' >> Results/quasar_replicate_results.txt
  done
fi

# Generate Noise Quality results
if [ ! -e Results/quasar_noise_quality_results.txt ]; then
  echo "Noise Quality"
  mkdir -p Noise/Results
  A=($(eval grep -v dm6 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    B=(${SAM//_/ })
    GEN=${B[1]}
    RE=${B[3]}
    if [ ! -e Data/${GEN}_${RE}_bin_corrections.hdf5 ]; then
      ${RUN} python Scripts/find_quasar_noise_model.py Data/\*${GEN}\*${RE}\*bin.hcp Data/${GEN}_${RE}_bin_corrections.hdf5 40000
    fi
    if [ -e Data/${GEN}_${RE}_bin_corrections.hdf5 -a ! -e Noise/${SAM}.quasar ]; then
      echo "noise quality ${SAM}"
      ${RUN} python Scripts/QuasarNoise.py Data/${SAM}_bin.hcp Data/${GEN}_${RE}_bin_corrections.hdf5 Noise/${SAM}.quasar -r 1000000,200000,40000 -d 0 -n 0.0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.35,0.5,0.75
    fi
    if [ ! -e Noise/Results/${SAM}.txt ]; then
      echo "noise quality score ${SAM}"
      python Scripts/find_quasar_quality_score.py Noise/${SAM}.quasar Noise/Results/${SAM}.txt
    fi
  done

  A=($(eval grep dm6 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    B=(${SAM//_/ })
    GEN=${B[1]}
    RE=${B[3]}
    if [ ! -e Data/${GEN}_${RE}_bin_corrections.hdf5 ]; then
      ${RUN} python Scripts/find_quasar_noise_model.py Data/\*${GEN}\*${RE}\*bin.hcp Data/${GEN}_${RE}_bin_corrections.hdf5 4000
    fi
    if [ -e Data/${GEN}_${RE}_bin_corrections.hdf5 -a ! -e Noise/${SAM}.quasar ]; then
      echo "noise quality ${SAM}"
      ${RUN} python Scripts/QuasarNoise.py Data/${SAM}_bin.hcp Data/${GEN}_${RE}_bin_corrections.hdf5 Noise/${SAM}.quasar -r 100000,20000,4000 -d 0 -n 0.0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.35,0.5,0.75
    fi
    if [ ! -e Noise/Results/${SAM}.txt ]; then
      echo "noise quality score ${SAM}"
      python Scripts/find_quasar_quality_score.py Noise/${SAM}.quasar Noise/Results/${SAM}.txt
    fi
  done

  A=($(eval cat Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  echo -e "Sample\tGenome\tNoise\tResolution\tScore" > Results/quasar_noise_quality_results.txt
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    B=(${SAM//_/ })
    GEN=${B[1]}
    if [ -e Noise/Results/${SAM}.txt ]; then
      tail -n +2 Noise/Results/${SAM}.txt | awk -v SAM="${SAM}" -v GEN="${GEN}" '{print SAM,GEN,$1,$2,$3}' OFS='\t' >> Results/quasar_noise_quality_results.txt
    fi
  done
fi

# Generate Noise Replicate results
if [ ! -e Results/quasar_noise_replicate_results.txt ]; then
  echo "Noise Replicate"
  A=($(eval grep Rep1 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  echo -e "Sample1\tSample2\tGenome\tNoise\tResolution\tScore" > Results/quasar_noise_replicate_results.txt
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]/_Rep1/}
    B=(${SAM//_/ })
    GEN=${B[1]}
    if [ -e Noise/${SAM}_Rep1.quasar -a -e Noise/${SAM}_Rep2.quasar ]; then
      if [ ! -e Noise/Results/${SAM}_Rep1_rep.txt ]; then
        echo "noise replicate ${SAM}_Rep1"
        python Scripts/find_quasar_replicate_noise_pair_score.py Noise/${SAM}_Rep1.quasar Noise/${SAM}_Rep2.quasar Noise/Results/${SAM}_Rep1_rep.txt
      fi
      tail -n +2 Noise/Results/${SAM}_Rep1_rep.txt | awk -v GEN="$GEN" '{print $1,$2,GEN,$3,$4,$5}' OFS='\t' >> Results/quasar_noise_replicate_results.txt
      if [ ! -e Noise/Results/${SAM}_Rep2_rep.txt ]; then
        echo "noise replicate ${SAM}_Rep2"
        python Scripts/find_quasar_replicate_noise_pair_score.py Noise/${SAM}_Rep2.quasar Noise/${SAM}_Rep1.quasar Noise/Results/${SAM}_Rep2_rep.txt
      fi
      tail -n +2 Noise/Results/${SAM}_Rep2_rep.txt | awk -v GEN="$GEN" '{print $1,$2,GEN,$3,$4,$5}' OFS='\t' >> Results/quasar_noise_replicate_results.txt
    fi
  done
fi

# Generate Coverage Quality results
if [ ! -e Results/quasar_coverage_quality_results.txt ]; then
  echo "Coverage Quality"
  mkdir -p Coverage/Results
  echo -e "Sample\tGenome\tCoverage\tResolution\tScore" > Results/quasar_coverage_quality_results.txt
  A=($(eval grep -v dm6 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    C=(${SAM//_/ })
    GEN=${C[1]}
    if [ ! -e Coverage/${SAM}.quasar ]; then
      echo "coverage quality ${SAM}"
      ${RUN} hifive quasar -p Data/${SAM}_bin.hcp Coverage/${SAM}.quasar -r 1000000,200000,40000 -d 40000000,20000000,10000000,5000000,2000000,1000000
    fi
    if [ ! -e Coverage/Results/${SAM}.txt -a -e Coverage/${SAM}.quasar ]; then
      echo "coverage quality score ${SAM}"
      python Scripts/find_quasar_quality_score.py Coverage/${SAM}.quasar Coverage/Results/${SAM}.txt
    fi
    if [ -e Coverage/Results/${SAM}.txt ]; then
      tail -n +2 Coverage/Results/${SAM}.txt | awk -v SAM="$SAM" -v GEN="$GEN" '{print SAM,GEN,$1,$2,$3}' OFS='\t' >> Results/quasar_coverage_quality_results.txt
    fi
  done

  A=($(eval grep dm6 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]}
    C=(${SAM//_/ })
    GEN=${C[1]}
    if [ ! -e Coverage/${SAM}.quasar ]; then
      echo "coverage quality ${SAM}"
      ${RUN} hifive quasar -p Data/${SAM}_bin.hcp Coverage/${SAM}.quasar -r 100000,20000,4000 -d 40000000,20000000,10000000,5000000,2000000,1000000,500000,250000,125000
    fi
    if [ ! -e Coverage/Results/${SAM}.txt -a -e Coverage/${SAM}.quasar ]; then
      echo "coverage quality score ${SAM}"
      python Scripts/find_quasar_quality_score.py Coverage/${SAM}.quasar Coverage/Results/${SAM}.txt
    fi
    if [ -e Coverage/Results/${SAM}.txt ]; then
      tail -n +2 Coverage/Results/${SAM}.txt | awk -v SAM="$SAM" -v GEN="$GEN" '{print SAM,GEN,$1,$2,$3}' OFS='\t' >> Results/quasar_coverage_quality_results.txt
    fi
  done
fi

# Generate Coverage Replicate results
if [ ! -e Results/quasar_coverage_replicate_results.txt ]; then
  echo "Coverage Replicate"
  echo -e "Sample1\tSample2\tGenome\tCoverage\tResolution\tScore" > Results/quasar_coverage_replicate_results.txt
  A=($(eval grep Rep1 Data/sample_key.txt | awk '{ print $1 }'))
  LENGTH=${#A[@]}
  for (( I=0; I<${LENGTH}; I++ )); do
    SAM=${A[${I}]/_Rep1/}
    B=(${SAM//_/ })
    GEN=${B[1]}
    if [ -e Coverage/${SAM}_Rep1.quasar -a -e Coverage/${SAM}_Rep2.quasar ]; then
      if [ ! -e Coverage/Results/${SAM}_rep.txt ]; then
        echo "coverage replicate ${SAM}"
        python ../../../Scripts/Python/find_quasar_replicate_pair_score.py Coverage/${SAM}_Rep1.quasar Coverage/${SAM}_Rep2.quasar Coverage/Results/${SAM}_rep.txt
      fi
      if [ -e Coverage/Results/${SAM}_rep.txt ]; then
        tail -n +2 Coverage/Results/${SAM}_rep.txt | awk -v GEN="${GEN}" '{print $1,$2,GEN,$3,$4,$5}' OFS='\t' >> Results/quasar_coverage_replicate_results.txt
      fi
    fi
  done
fi

# Generate heterogeneity results
if [ ! -e Results/quasar_heterogeneity_results.txt ]; then
  echo "Quasar Heterogeneity"
  mkdir -p Heterogeneity/Results
  A=($(eval cat Data/heterogeneity_pairs.txt))
  LENGTH=${#A[@]}
  echo -e "Sample1\tSample2\tBalance\tCoverage\tResolution\tScore" > Results/quasar_heterogeneity_results.txt
  for (( I=0; I<${LENGTH}; I+=2 )); do
    SAM1=${A[${I}]}
    SAM2=${A[${I}+1]}
    if [ -e Coverage/Results/${SAM1}.txt ]; then
      grep 20000000 Coverage/Results/${SAM1}.txt | grep -v -P '10000\t' | awk -v SAM1=${SAM1} -v SAM2=${SAM2} -v BAL='1.0' '{print SAM1,SAM2,BAL,$1,$2,$3}' OFS='\t' >> Results/quasar_heterogeneity_results.txt
    fi
    if [ -e Coverage/Results/${SAM2}.txt ]; then
      grep 20000000 Coverage/Results/${SAM2}.txt | grep -v -P '10000\t' | awk -v SAM1=${SAM1} -v SAM2=${SAM2} -v BAL='0.0' '{print SAM1,SAM2,BAL,$1,$2,$3}' OFS='\t' >> Results/quasar_heterogeneity_results.txt
    fi
    for J in 0.2 0.4 0.6 0.8; do
      if [ ! -e Heterogeneity/${SAM1}_${SAM2}_${J}.quasar ]; then
        echo "heterogeneity ${SAM1} ${SAM2} ${J}"
        ${RUN} python Scripts/QuasarPseudo.py Data/${SAM1}_bin.hcp Data/${SAM2}_bin.hcp Heterogeneity/${SAM1}_${SAM2}_${J}.quasar -r 1000000,200000,40000 -d 20000000 -b ${J}
      fi
      if [ ! -e Heterogeneity/Results/${SAM1}_${SAM2}_${J}.txt ]; then
        echo "heterogeneity score ${SAM1} ${SAM2} ${J}"
        python Scripts/find_quasar_quality_score.py Heterogeneity/${SAM1}_${SAM2}_${J}.quasar Heterogeneity/Results/${SAM1}_${SAM2}_${J}.txt
      fi
      if [ -e Heterogeneity/Results/${SAM1}_${SAM2}_${J}.txt ]; then
        tail -n +2 Heterogeneity/Results/${SAM1}_${SAM2}_${J}.txt | awk -v SAM1=${SAM1} -v SAM2=${SAM2} -v BAL=${J} '{print SAM1,SAM2,BAL,$1,$2,$3}' OFS='\t' >> Results/quasar_heterogeneity_results.txt
      fi
    done
  done
fi

# Generate Figure 1
mkdir -p Figures
for I in Encode_hg38_A549_HindIII_Rep1 GSE35156_mm10_ES_HindIII_Rep1 GSE60494_mm10_ES_MboI_Rep2; do
  rm -f Results/${I}.quasar
  if [ -e Quasar/${I}.quasar ]; then
    ln -s ${PWD}/Quasar/${I}.quasar Results/
  else
    ln -s ${PWD}/Results/${I}_subset.quasar Results/${I}.quasar
  fi
done
python Scripts/Figure_1.py

# Generate Figure 2
python Scripts/Figure_2.py

# Generate Figure 3
python Scripts/Figure_3.py

# Generate Figure 4
python Scripts/Figure_4.py

# Generate Figure 5
I=(GSM1267201 GSM1267203 Encode_hg38_LNCaP-clone-FGC_HindIII_Rep2)
J=(GSE52457_hg38_H1-MesenchymalStem_HindIII_Rep2 GSE52457_hg38_H1-NPC_HindIII_Rep2 Encode_hg38_LNCaP-clone-FGC_HindIII_Rep2)
for (( K=0; K<${#I[*]}; K++ )); do
  SAM=${I[${K}]}
  ALIAS=${J[${K}]}
  C=(${ALIAS//_/ })
  if [ ! -e Data/${SAM}.hcd ]; then
    wget --no-check-certificate -O Data/${SAM}.hcd  https://bx.bio.jhu.edu/data/quasar/${SAM}.hcd
    wget --no-check-certificate -O Data/${SAM}_bin.hcp  https://bx.bio.jhu.edu/data/quasar/${SAM}_bin.hcp
      if [ ${SAM} != ${ALIAS} ]; then
        ln -s ${PWD}/Data/${SAM}_bin.hcp Data/${ALIAS}_bin.hcp
      fi
  fi
  if [ ! -e Data/${C[1]}_${C[3]}.fends ]; then
    wget --no-check-certificate -O Data/${C[1]}_${C[3]}.fends  https://bx.bio.jhu.edu/data/quasar/${C[1]}_${C[3]}.fends
  fi
done
python Scripts/Figure_5.py

# Generate Figure S1
python Scripts/Figure_S1.py

# Generate Figure S2
python Scripts/Figure_S2.py

# Generate Figure S3
python Scripts/Figure_S3.py

# Generate Figure S4
python Scripts/Figure_S4.py

# Generate Table S1
python Scripts/Table_S1.py
