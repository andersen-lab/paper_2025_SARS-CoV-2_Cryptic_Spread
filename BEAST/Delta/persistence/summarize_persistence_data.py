import argparse
import pandas as pd
import baltic as bt
import glob
import os
import numpy as np
import arviz as az

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--location', type=str, required=False, help='Location of Interest', default="Iraq")
	parser.add_argument('--lineage', type=str, required=False, help='Lineage of Interest', default="delta5")
	parser.add_argument('--tree', type=str, required=False, help='Path to the tree file', default="/raid/gp/KUR/Datasets/Delta/dataset-3/xml3/dup_run/geo/div_EU_AS_2/tree.nexus")
	return parser.parse_args()


def summarise_group(group):
	anc_eval_diff = (group['ancestralTime'] - group['evaluationTime']).iloc[0]
	# Number of persistent lineages at eval time
	persistentLineagesAtEval = (group['persistenceTime'] > anc_eval_diff).sum()
	# Number of lineages from new introductions at eval time
	introductionLineagesAtEval = (group['persistenceTime'] <= anc_eval_diff).sum()
	# independenceTime used to filter to "unique" lineages
	persistentsFromUnique = ((group['independenceTime'] > anc_eval_diff) & (group['persistenceTime'] > anc_eval_diff)).sum()
	introductionsFromUnique = ((group['independenceTime'] > anc_eval_diff) & (group['persistenceTime'] <= anc_eval_diff)).sum()
	total_unique = persistentsFromUnique + introductionsFromUnique
	propPersistentFromUnique = persistentsFromUnique / total_unique if total_unique > 0 else 0
	return pd.Series({
        'anc_eval_diff': anc_eval_diff,
        'persistentLineagesAtEval': persistentLineagesAtEval,
        'introductionLineagesAtEval': introductionLineagesAtEval,
        'persistentsFromUnique': persistentsFromUnique,
        'introductionsFromUnique': introductionsFromUnique,
        'propPersistentFromUnique': propPersistentFromUnique,
		'total_unique': total_unique
    })

def summarise_hpd_lower(x):
	return az.hdi(np.array(x), hdi_prob=0.95)[0]

def summarise_hpd_upper(x):
    return az.hdi(np.array(x), hdi_prob=0.95)[1]

def main():
	args = parse_args()

	# BAN = "2021-04-27"
     
	location = args.location
	lineage = args.lineage
	tree = args.tree
	# location = "Iraq"
	# lineage = "delta5"
	# tree = "/raid/gp/KUR/Datasets/Delta/dataset-3/xml3/dup_run/geo/div_EU_AS_2/tree.nexus"
	
	print(f"Location: {location}")
	print(f"Lineage: {lineage}")

	out_folder = f"./output/{lineage}/"
	out_files = glob.glob(os.path.join(out_folder, "**/*.tsv"), recursive=True)

	print(f"Output files found are {out_files}")

    # Read and concatenate data
	df = pd.read_csv(out_files[0], engine='pyarrow')
	
	# Remove non-ASCII characters
	df['stateAtEvaluationTime'] = df['stateAtEvaluationTime'].str.encode('ascii', 'ignore').str.decode('ascii')

	# mrsd = bt.decimalDate('2022-01-30')
	myTree = bt.loadNexus(tree)

	# mrsd: most recent sampling date
	mrsd = max([i.absoluteTime for i in myTree.Objects])
	print(mrsd)
	print(bt.calendarDate(mrsd))

	# filter to include only Iraq
	df_filtered = df[df['stateAtEvaluationTime'] == location]

	# summarize data
	grouped = df_filtered.groupby(['treeId', 'evaluationTime', 'ancestralTime', 'stateAtEvaluationTime'])
	df_prop = grouped.apply(summarise_group).reset_index()

	# Convert times to dates
	df_prop['evaluationTime'] = mrsd - df_prop['evaluationTime']
	df_prop['ancestralTime'] = mrsd - df_prop['ancestralTime']
	df_prop['mean_time'] = (df_prop['ancestralTime'] + df_prop['evaluationTime']) / 2

	print(df_prop['evaluationTime'].value_counts())
	print(df_prop['ancestralTime'].value_counts())

	df_prop['treeId'] = df_prop['treeId'].astype('category')

	summary_stats = df_prop.groupby(['ancestralTime', 'evaluationTime', 'mean_time']).agg(
        # Proportion of persistent lineages
    	propPersistentFromUnique_lower=('propPersistentFromUnique', summarise_hpd_lower),
    	propPersistentFromUnique_upper=('propPersistentFromUnique', summarise_hpd_upper),
    	propPersistentFromUnique_median=('propPersistentFromUnique', np.median),
        # Number of persistent lineages
    	persistentsFromUnique_lower=('persistentsFromUnique', summarise_hpd_lower),
    	persistentsFromUnique_upper=('persistentsFromUnique', summarise_hpd_upper),
    	persistentsFromUnique_median=('persistentsFromUnique', np.median),
        # Number of unique introductions
    	introductionsFromUnique_lower=('introductionsFromUnique', summarise_hpd_lower),
    	introductionsFromUnique_upper=('introductionsFromUnique', summarise_hpd_upper),
    	introductionsFromUnique_median=('introductionsFromUnique', np.median),
	).reset_index()
     
	summary_stats.to_csv(f"./output/{lineage}/{location}_summary_stats.tsv", sep='\t', index=False)
    
	# %%
	# Proportion of lineages at evalTime
	df_prop['propDescendantsFromIntroductionsAtEval'] = (
        df_prop['introductionLineagesAtEval'] /
        (df_prop['introductionLineagesAtEval'] + df_prop['persistentLineagesAtEval'])
    )

	# print(df_prop['propDescendantsFromIntroductionsAtEval'])

	descendants_summary = df_prop.groupby(['ancestralTime', 'evaluationTime', 'mean_time']).agg(
		propDescendantsFromIntroductionsAtEval_lower=('propDescendantsFromIntroductionsAtEval', summarise_hpd_lower),
		propDescendantsFromIntroductionsAtEval_upper=('propDescendantsFromIntroductionsAtEval', summarise_hpd_upper),
		propDescendantsFromIntroductionsAtEval_median=('propDescendantsFromIntroductionsAtEval', np.median)
    )
	descendants_summary = descendants_summary.apply(pd.Series)
	descendants_summary = descendants_summary.reset_index()
     
	descendants_summary.to_csv(f"./output/{lineage}/{location}_descendants_summary.tsv", sep='\t', index=False)
     
	# %%
	# introductionLineagesAtEval
	introductionLineagesAtEvalStats = df_prop.groupby(['ancestralTime', 'evaluationTime', 'mean_time']).agg(
		introductionLineagesAtEval_lower=('introductionLineagesAtEval', summarise_hpd_lower),
		introductionLineagesAtEval_upper=('introductionLineagesAtEval', summarise_hpd_upper),
		introductionLineagesAtEval_median=('introductionLineagesAtEval', np.median)
	)

	introductionLineagesAtEvalStats = introductionLineagesAtEvalStats.apply(pd.Series)
     
	introductionLineagesAtEvalStats = introductionLineagesAtEvalStats.reset_index()

	introductionLineagesAtEvalStats.to_csv(f"./output/{lineage}/{location}_introductionLineagesAtEvalStats.tsv", sep='\t', index=False)
# %%
if __name__ == '__main__':
	main()
