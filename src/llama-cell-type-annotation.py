import os
from typing import Union, List, Dict
import pandas as pd
from hugchat import hugchat
from hugchat.login import Login

def huggingchatcelltype(
    input: Union[pd.DataFrame, Dict[str, List[str]]],
    model: int = 0,
    tissuename: str = None,
    topgenenumber: int = 20,
    add_info: str = None,
    username: str = None,
    password: str = None
) -> Union[str, List[str]]:
    """
    Annotate cell types using Llama 3 model via HuggingChat.
    
    Parameters:
    -----------
    input : Union[pd.DataFrame, Dict[str, List[str]]]
        Either a Seurat-style differential gene DataFrame or a dictionary of gene lists
    model : int
        Index of the model to use.
    tissuename : str, optional
        Name of the tissue being analyzed
    topgenenumber : int, default 10
        Number of top differential genes to use if input is a DataFrame
    add_info : str, optional
        Additional context to help with cell type identification
    username : str, optional
        HuggingChat username (if not provided, will look for environment variable)
    password : str, optional
        HuggingChat password (if not provided, will look for environment variable)
    
    Returns:
    --------
    Union[str, List[str]]
        If no login credentials are provided, returns the prompt
        Otherwise, returns a list of cell type annotations
    """
    # Intro message
    print("Welome to hugchat cell-type annotation tool")

    # Check for login credentials
    if username is None:
        username = os.getenv('HUGGINGCHAT_USERNAME')
    if password is None:
        password = os.getenv('HUGGINGCHAT_PASSWORD')
    
    # Process input 
    if isinstance(input, dict):
        # If input is a dictionary, convert gene lists to comma-separated strings
        processed_input = {k: ','.join(v) for k, v in input.items()}
    elif isinstance(input, pd.DataFrame):
        # If input is a DataFrame, filter significant genes and get top genes per cluster
        input = input[input['avg_log2FC'] > 0]
        processed_input = {}
        for cluster in input['cluster'].unique():
            cluster_genes = input[input['cluster'] == cluster]['gene']
            processed_input[cluster] = ','.join(cluster_genes[:topgenenumber])
    else:
        raise ValueError("Input must be a pandas DataFrame or a dictionary of gene lists")
    
    # Construct prompt
    prompt = (
        f"Identify cell types of {tissuename or 'unknown'} cells using the following markers separately for each row. "
        "Only provide the cell type name. Do not show numbers before the name. "
        "Some could be a mixture of multiple cell types. "
        f"{add_info or ''}\n" +
        "\n".join([f"{cluster}: {genes}" for cluster, genes in processed_input.items()])
    )
    
    # If no login credentials, return prompt
    if not username or not password:
        print("Note: HuggingChat login credentials not found: returning the prompt itself.")
        return prompt
    
    # Login to HuggingChat
    try:
        sign = Login(username, password)
        cookies = sign.login()
        
        # Create ChatBot 
        chatbot = hugchat.ChatBot(cookies=cookies)

        # Select model - this can be done by selecting the index of the model from the list of models
        chatbot.switch_llm(model)

        # Print list of available models
        available = chatbot.get_remote_llms()
        print(f'The selected model is : {available[model].name}')
        print('Note that other models are available and can be selected by changing the model parameter :\n')
        for i, model in enumerate(available):
            print(f'{i} : {model.name}')
        
        
        
        # Send prompt and get response
        cell_type_annotations = []
        
        # Process in batches to handle large numbers of clusters
        batch_size = 100
        clusters = list(processed_input.keys())
        
        for i in range(0, len(clusters), batch_size):
            batch_clusters = clusters[i:i+batch_size]
            batch_prompt = (
                f"Identify cell types of {tissuename or 'unknown'} cells using the following markers separately for each row. "
                "Only provide the cell type name. Do not show numbers before the name. "
                f"{add_info or ''}\n" +
                "\n".join([f"{processed_input[cluster]}" for cluster in batch_clusters])
            )
            
            response = chatbot.chat(batch_prompt).wait_until_done()
            batch_annotations = response.split('\n')
            
            # Trim and validate annotations
            batch_annotations = [
                annotation.strip() 
                for annotation in batch_annotations 
                if annotation.strip() and not annotation.startswith(tuple('0123456789'))
            ]
            
            # Ensure we have one annotation per cluster in the batch
            if len(batch_annotations) != len(batch_clusters):
                print(f"Warning: Unexpected number of annotations in batch. Expected {len(batch_clusters)}, got {len(batch_annotations)}")
                # Pad or truncate to match
                batch_annotations = batch_annotations[:len(batch_clusters)] + ['Unknown'] * max(0, len(batch_clusters) - len(batch_annotations))
            
            cell_type_annotations.extend(batch_annotations)
        
        print('Annotation done !')
        print('Note: It is always recommended to check the results returned by this function in case of AI hallucination, before proceeding with downstream analysis.')
        
        return cell_type_annotations
    
    except Exception as e:
        print(f"Error in HuggingChat interaction: {e}")
        return prompt

# Example usage
if __name__ == "__main__":
    # Example with dictionary input
    gene_dict = {
        'cluster1': ['CD4', 'CD3D'],
        'cluster2': ['CD14']
    }
    
    # Set HUGGINGCHAT_USERNAME and HUGGINGCHAT_PASSWORD environment variables before running
    annotations = huggingchatcelltype(
        input=gene_dict, 
        model=1,
        tissuename='human PBMC', 
        add_info=''
    )
    print(annotations)
