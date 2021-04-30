#! /usr/bin/env awk

NR==FNR{
    tpm[$1"\t"$2]=$3
}
NR>FNR && FNR>1{
    split($7, inclusion, ",");
    split($8, exclusion, ",");
    inclusion_sample1=0
    inclusion_sample2=0
    inclusion_sample3=0
    inclusion_sample4=0
    exclusion_sample1=0
    exclusion_sample2=0
    exclusion_sample3=0
    exclusion_sample4=0
    sample1="C08leaf_0h"
    sample2="C08leaf_1h"
    sample3="C08root_0h"
    sample4="C08root_1h"
    for(i=1;i<=length(inclusion);i++){
        inclusion_sample1 += tpm[inclusion[i]"\t"sample1]
        inclusion_sample2 += tpm[inclusion[i]"\t"sample2]
        inclusion_sample3 += tpm[inclusion[i]"\t"sample3]
        inclusion_sample4 += tpm[inclusion[i]"\t"sample4]
    }
    for(i=1;i<=length(exclusion);i++){
        exclusion_sample1 += tpm[exclusion[i]"\t"sample1]
        exclusion_sample2 += tpm[exclusion[i]"\t"sample2]
        exclusion_sample3 += tpm[exclusion[i]"\t"sample3]
        exclusion_sample4 += tpm[exclusion[i]"\t"sample4]
    }
    ## if(inclusion_sample1==0 && exclusion_sample1==0){
    ##     psi_sample1="NA"
    ## }else{
    ##     psi_sample1=inclusion_sample1/(inclusion_sample1+exclusion_sample1)
    ## }
    ## if(inclusion_sample2==0 && exclusion_sample2==0){
    ##     psi_sample2="NA"
    ## }else{
    ##     psi_sample2=inclusion_sample2/(inclusion_sample2+exclusion_sample2)
    ## }
    ## if(inclusion_sample3==0 && exclusion_sample3==0){
    ##     psi_sample3="NA"
    ## }else{
    ##     psi_sample3=inclusion_sample3/(inclusion_sample3+exclusion_sample3)
    ## }
    ## if(inclusion_sample4==0 && exclusion_sample4==0){
    ##     psi_sample4="NA"
    ## }else{
    ##     psi_sample4=inclusion_sample4/(inclusion_sample4+exclusion_sample4)
    ## }
    inclusion_expr = inclusion_sample1"|"inclusion_sample2"|"inclusion_sample3"|"inclusion_sample4
    exclusion_expr = exclusion_sample1"|"exclusion_sample2"|"exclusion_sample3"|"exclusion_sample4
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"inclusion_expr"\t"exclusion_expr
}
