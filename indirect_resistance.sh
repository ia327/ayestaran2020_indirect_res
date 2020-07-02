#!/bin/bash
while [[ $# -gt 0 ]]
  do key="$1"
     case $key in
        -h)  printf "############   HELP   ###############\nOPTIONS\n"
            printf "\t-h\tThis Help\n"
            printf "\t-C\tForce the generation of combined_data.RData\n"
            printf "\t-A\tRun tissue-specific ANOVA\n"
            printf "\t-g|--gdsc\tUse GDSC dataset\n"
            printf "\t-c|--ctrp\tUse CTRP dataset\n"
            exit 0
            ;;
        -C)
        force_combined_data=1
        shift
        ;;
        -A)
        run_ANOVA=1
        shift
        ;;
        -g|--gdsc)
        dataset='gdsc'
        shift
        ;;
        -c|--ctrp)
        dataset='ctrp'
        shift
        ;;
    esac
  done


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

upperdataset=`echo $dataset | tr '[:lower:]' '[:upper:]'`

if [ "$dataset" != "gdsc" ] && [ "$dataset" != "ctrp" ]
    then
        printf "Error: choose a valid dataset: GDSC or CTRP\n"
        exit 0
    else
        printf "Analysing data from $upperdataset data set\n"
fi

# Run from main directory
# if [ ! -f "$DIR/data/combined_data.RData" ]
#     then
#         printf "Combined data file is missing, a new one will be created\n"
#         printf "Next time, use -C to create the combined data file again\n"
#         force_combined_data=1
# fi

if [ "$force_combined_data" == 1 ]
    then
        printf "Creating combined data file...\n"
        Rscript --vanilla "$DIR/scripts/create_combined_data.R" "$DIR"
fi

mkdir -p "$DIR/results"
mkdir -p "$DIR/results/$dataset/"
mkdir -p "$DIR/results/$dataset/ANOVA_out/"

if [ "$run_ANOVA" == 0 ] && [ ! -d "$DIR/results/ANOVA_out" ]
    then printf "Error: results from ANOVA are missing, try using -A\n"
    exit 0
fi

if [ "$run_ANOVA" == 1 ]
    then
        printf "Running tissue specific ANOVA...\n"
        Rscript --vanilla "$DIR/scripts/run_all_anovas.R" "$DIR" "$dataset" \
        "$DIR/results/$dataset/ANOVA_out/"
fi

mkdir -p "$DIR/plots"
mkdir -p "$DIR/plots/$dataset"


printf "Filtering ANOVA output...\n"
Rscript --vanilla "$DIR/scripts/anova_out_filter.R" "$DIR" "$dataset" \
"$DIR/results/$dataset/ANOVA_out/" "$DIR/results/$dataset/"

printf "Running variance change analysis...\n"
Rscript --vanilla "$DIR/scripts/variance_analysis.R" "$DIR" "$dataset" \
"$DIR/results/$dataset/" "$DIR/results/$dataset/"

printf "Looking for putative markers within found outliers...\n"
Rscript --vanilla "$DIR/scripts/seq_cn_lookup.R" "$DIR" "$dataset" \
"$DIR/results/$dataset/" "$DIR/results/$dataset/"

printf "Done\n"
