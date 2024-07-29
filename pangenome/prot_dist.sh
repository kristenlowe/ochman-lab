# number of proteins in each strain
for i in {1..72}; do
    grep ">" ECOR_${i}_protein.faa | wc -l
done

4339
4900
4544
4198
4781
4126
grep: ECOR_7_protein.faa: No such file or directory
0
4469
4870
4321
4823
4791
4225
4574
4555
4205
4088
4210
4123
4169
4139
4137
grep: ECOR_23_protein.faa: No such file or directory
0
4773
4397
4343
4567
4578
4593
4495
4849
4418
4420
4563
4741
4909
5292
4845
4912
4823
4869
4833
4896
4941
4381
4840
4578
4927
4888
5225
4779
4637
4684
4600
4686
4536
4854
4541
4391
4631
4560
4812
4717
4715
4497
4340
4397
4658
4260
4733
4527
4375

# get list of all protein names in ECOR_1_protein.faa
grep ">" ECOR_1_protein.faa | sed 's/^>//; s/ .*$//'


# protein ECOR_01_A_WP_001308243.1 is present in these strains
awk -F '\t' '($3>70&&$5>70&&$16<0.001)' ECOR_1_vs_ECOR.tsv | grep "ECOR_01_A_WP_001308243.1" | cut -f2 | cut -f1-3 -d "_" | sort -u

