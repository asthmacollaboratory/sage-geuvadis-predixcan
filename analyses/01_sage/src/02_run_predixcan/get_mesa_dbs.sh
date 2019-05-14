#!/usr/bin/env bash

db_sfx="imputed_10_peer_3_pcs_v2.db"

# get DBs for AFA, AFHI, CAU
for pop in AFA AFHI CAU; do
    echo "downloading DB for $pop"
    db="${pop}_${db_sfx}"
    wget "https://s3.amazonaws.com/predictdb2/contributed/MESA-2018-05-v2/${db}"
done

# do ALL manually since it has a typo in file path
wget https://s3.amazonaws.com/predictdb2/contributed/MESA-2018-05-v2/ALL_imputed_10_peer_3_pcs__v2.db
mv ALL_imputed_10_peer_3_pcs__v2.db ALL_${db_sfx}

for pop in AFA AFHI CAU ALL; do
    echo "parsing DB file for $pop"
    db="${pop}_${db_sfx}"
    outfile="mesa.${pop}.test.metrics.txt"
    sqlite3 -header -separator "," "./${db}" "select gene, genename, test_R2_avg, rho_avg from extra;" > $outfile 
done
