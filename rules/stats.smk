

rule agg_all_stats:
    input: expand_rows_with_lookup(paths.basecalls.agg_stats, basecalls_df, is_aggreated=True) 
    output: paths.stats.agg_stat
    # conda: 
    run:
        import pandas as pd 
        from pathlib import Path

        def flowcell(x):
            return str(Path(x).name).split('.')[0] 

        df = pd.concat([pd.read_csv(f).assign(flowcell=flowcell(f)) for f in input])
        df.to_csv(output[0])
