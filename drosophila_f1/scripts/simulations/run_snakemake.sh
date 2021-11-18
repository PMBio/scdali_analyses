snakemake --rerun-incomplete --latency-wait 30 --restart-times 3 --cluster "bsub -n {threads} -R rusage[mem={resources.mem_mb}] -e lsf_reports/T-%J.err -o lsf_reports/T-%J.out -q medium"  -j 30 #-T 2
