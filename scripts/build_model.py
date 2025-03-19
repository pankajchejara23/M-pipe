import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.metrics import accuracy_score, precision_recall_curve, auc

import matplotlib.pyplot as plt
from sklearn.linear_model import LassoCV
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split



def get_specific_label(l, t):
    """
    This function parse the taxonomic assignment lable and fetch the specified information (e.g., kingdom, family)

    Args:
    -----
        l (str): string of taxonomy
        t (str): string specifying the requested information (e.g., kingdom, family, genus, etc.)

    Returns:
    -----
        str: requested heirarichal info 
    """
    taxa_order = {'kingdom':0,'phylum':1,'class':2,'order':3,'family':4,'genus':5,'species':6}

    try:
        specific_label = l.split(';')[taxa_order[t]]

        return specific_label.strip().split('__')[1]
    except:
        return 'Unknown'


def getTaxa(taxa_df, feature, t):
    """
    This function obtains taxonomy of given taxon.

    Args:
    ----
        taxa (DataFrame) - Pandas dataframe containing taxonomy of features
        feature (str) - feature id 
        t (str): string specifying the requested information (e.g., kingdom, family, genus, etc.)

    Returns:
    ----
        str: Taxonomy associated with feature at level t
    """

    row = taxa_df.loc[taxa['Feature ID'] == feature,'Taxon'].to_list()
    if len(row) != 0:
        taxonomy = row[0]
        return get_specific_label(taxonomy,t)
    else:
        print('Feature not exists')
        return None

def aggregate(otu_df, rank):
    """
    This function aggregate counts at specified rank level.

    Args:
    ----
        otu_df (DataFrame) - Pandas Dataframe containing microbial abundance count data.
        rank (str) - Rank level at which to perform aggregation

    Returns:
    ----
        DataFrame - Aggregated dataframe
    """
    # Extract rank info for each feature
    otu_df[rank] = otu_df['#OTU ID'].apply(lambda x: getTaxa(taxa,x,rank))

    # Drop OTU ID
    otu = otu_df.drop(['#OTU ID'],axis=1)

    # Groupby on Rank
    otu_df_sum = otu.groupby([rank]).sum()
    otu_df_sum.reset_index(inplace=True)
    rank_labels = otu_df_sum[rank].to_list()
    otu_df_sum.drop([rank],axis=1,inplace=True)

    # Transpose and set rank column as header
    otu_df_T = otu_df_sum.T
    otu_df_T.columns = rank_labels
    
    return np.log1p(otu_df_T)

def attach_labels(otu,meta,index_col, target_col, target_group):
    """
    This function attach target labels (e.g., disease type) to otu data.

    Args:
    ----
        otu (DataFrame) - Pandas dataframe containing abundance count
        meta (DataFrame) - Metadata of samples
        index_col (str) - Name of column containing sample id
        target_col (str) - Name of target column

    Returns:
    ----
        DataFrame - A pandas dataframe containing otu and target labels
    
    """
    meta_ind = meta.set_index([index_col])[target_col]
    data = pd.concat([otu,meta_ind],axis=1)
    data['class'] = data[target_col].apply(lambda x: 1 if x == target_group else 0)
    return data.drop([target_col],axis=1)

def get_data_from_pipeline_output(otu, meta, rank, index_col, target_col, target_group):
    """
    This function takes output files from upstream pipeline and processes them to prepare for ML step.

    Args:
    ----
        otu (DataFrame) - Pandas dataframe containing abundance count
        meta (DataFrame) - Metadata of samples
        rank (str) - Rank level at which to perform aggregation
        index_col (str) - Name of column containing sample id
        target_col (str) - Name of target column

    Returns:
    ----
        DataFrame - Dataset
    
    """

    otu_agg = aggregate(otu,rank)
    data = attach_labels(otu_agg,meta,index_col,target_col,target_group)
    return data.drop(['class'],axis=1), data['class']

def lasso_cv_auc_plot(X, y, n_alphas=100, cv=5, random_state=42):
    """
    Perform LassoCV on X and y, compute AUC, and plot the ROC curve.

    Parameters:
    - X: Feature matrix (n_samples, n_features).
    - y: Target labels (n_samples,).
    - n_alphas: Number of alpha values to try in LassoCV.
    - cv: Number of cross-validation folds.
    - random_state: Random seed for reproducibility.

    Returns:
    - auc: AUC score.
    - Plot: ROC curve.
    """
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=random_state)

    # Perform LassoCV
    lasso_cv = LassoCV(n_alphas=n_alphas, max_iter=10000, cv=cv, random_state=random_state)
    lasso_cv.fit(X_train, y_train)

    # Predict probabilities (for AUC calculation)
    y_pred = lasso_cv.predict(X_test)

    # Compute AUC
    auc = roc_auc_score(y_test, y_pred)
    print(f"AUC: {auc:.2f}")

    # Plot ROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f"LassoCV (AUC = {auc:.2f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Random model")
    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")
    plt.title("ROC Curve for LassoCV")
    plt.legend()
    plt.grid()
    plt.show()

    return auc


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Microbiome Machine Learning Script")
    parser.add_argument("otu", type=str, help="Path to OTU abundance TSV file")
    parser.add_argument("taxa", type=str, help="Path to taxonomy TSV file")
    parser.add_argument("meta", type=str, help="Path to metadata TSV file")
    parser.add_argument("rank", type=str, help="Taxonomy rank to aggregate by (e.g., 'Genus')")
    parser.add_argument("index_col", type=str, help="Column name in metadata containing sample identifiers")
    parser.add_argument("target_col", type=str, help="Column name in metadata containing target labels")
    parser.add_argument("target_group", type=str, help="Group in target_col to treat as class 1")
    parser.add_argument("output", type=str, help="Output file to save auc plot")

    # Parse arguments
    args = parser.parse_args()

    # Load data
    otu = pd.read_csv(args.otu, sep='\t',skiprows=1)
    taxa = pd.read_csv(args.taxa, sep='\t')
    meta = pd.read_csv(args.meta,sep='\t')

    # Process data for modeling
    X,y = get_data_from_pipeline_output(otu,meta,args.rank,args.index_col,args.target_col,args.target_group)

    # Training LASSO and plotting AUC curve
    lasso_cv_auc_plot(X,y,plot_file)