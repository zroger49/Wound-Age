#!bin/bash

######## Run ORSUM (Directional)#############

# Control vs Early
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Up_Wound_Early --files "results_new/03_enrichement/Control - Early/orsum.go.bp.up.results.csv" --outputFolder "results_new/03_enrichement/Control - Early/ORSUM_bp_up"

# Control vs Intermediate
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Up_Wound_Intermediary --files "results_new/03_enrichement/Control - Intermediary/orsum.go.bp.up.results.csv" --outputFolder "results_new/03_enrichement/Control - Intermediary/ORSUM_bp_up"
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Down_Wound_Intemediary --files "results_new/03_enrichement/Control - Intermediary/orsum.go.bp.down.results.csv" --outputFolder "results_new/03_enrichement/Control - Intermediary/ORSUM_bp_down"

# Control vs Late
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Up_Wound_Late --files "results_new/03_enrichement/Control - Late/orsum.go.bp.up.results.csv" --outputFolder "results_new/03_enrichement/Control - Late/ORSUM_bp_up"
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Down_Wound_Late --files "results_new/03_enrichement/Control - Late/orsum.go.bp.down.results.csv" --outputFolder "results_new/03_enrichement/Control - Late/ORSUM_bp_down"

# Early vs Intemediary
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Down_Wound_Early_Intemediary --files "results_new/03_enrichement/Early vs Intermediary/orsum.go.bp.down.results.csv" --outputFolder "results_new/03_enrichement/Early vs Intermediary/ORSUM_bp_down"

# Intermediary vs Late
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Up_Wound_Intermediary_Late --files "results_new/03_enrichement/Intermediary vs Late/orsum.go.bp.up.results.csv" --outputFolder "results_new/03_enrichement/Intermediary vs Late/ORSUM_bp_up"
python3 code/orsum/orsum.py --gmt data/hsapiens.GO.BP.name.gmt --fileAliases Down_Wound_Intermediary_Late --files "results_new/03_enrichement/Intermediary vs Late/orsum.go.bp.down.results.csv" --outputFolder "results_new/03_enrichement/Intermediary vs Late/ORSUM_bp_down"

# Join all ORSUM results into a single summary file

python3 code/orsum/orsum.py \
  --gmt data/hsapiens.GO.BP.name.gmt \
  --fileAliases Up_Early Up_Intermediary Down_Intermediary Up_Late Down_Late Down_Early_Intermediary Up_Intermediary_Late Down_Intermediary_Late \
  --files "results_new/03_enrichement/Control - Early/orsum.go.bp.up.results.csv" \
           "results_new/03_enrichement/Control - Intermediary/orsum.go.bp.up.results.csv" \
           "results_new/03_enrichement/Control - Intermediary/orsum.go.bp.down.results.csv" \
            "results_new/03_enrichement/Control - Late/orsum.go.bp.up.results.csv" \
            "results_new/03_enrichement/Control - Late/orsum.go.bp.down.results.csv" \
            "results_new/03_enrichement/Early vs Intermediary/orsum.go.bp.down.results.csv" \
            "results_new/03_enrichement/Intermediary vs Late/orsum.go.bp.up.results.csv" \
            "results_new/03_enrichement/Intermediary vs Late/orsum.go.bp.down.results.csv" \
  --outputFolder results_new/03_enrichement/Summary_All_Contrasts \
  --minTermSize 5 
