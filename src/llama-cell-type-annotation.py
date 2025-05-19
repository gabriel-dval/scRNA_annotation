import os
from typing import Union, List, Dict
import pandas as pd
from hugchat import hugchat
from hugchat.login import Login


# Util functions to load data from seurat output

def seuratmarker_to_dict(csv_file: str, 
                         topgenenumber: int=20,
                         logfc_threshold: float=0,
    ) -> Dict[str, List[str]]:
    '''Function to convert seurat FindMarkers output (in csv format) to a dictionary of gene lists
    per cluster.

    Args
    ----
    csv_file : str
        Path to the csv file containing the seurat FindMarkers output.
    topgenenumber : int
        Number of top genes to consider per cluster.
    logfc_threshold : float
        Log Fold-Change threshold for significance.
    

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene lists per cluster.
    '''
    # Read the csv file
    input = pd.read_csv(csv_file)
    print(input)
    
    # Filter significant genes and most differentially expressed genes
    input = input[input['avg_log2FC'] > logfc_threshold]
    processed_input = {}
    for cluster in input['cluster'].unique():
            cluster_genes = input[input['cluster'] == cluster]['gene']
            processed_input[cluster] = ','.join(cluster_genes[:topgenenumber])

    return processed_input
    

def list_to_csv(data: List[str], output_file: str) -> None:
    """
    Save a list of elements into a single-column CSV file with the column name 'Annotation'.

    Args:
    ----
    data : List[str]
        List of elements to save.
    output_file : str
        Path to the output CSV file.

    Returns:
    -------
    None
    """
    df = pd.DataFrame(data, columns=['Annotation'])
    df.to_csv(output_file, index=False)


# Main function

def huggingchatcelltype(
    input: Union[str, Dict[str, List[str]]],
    output_dir: str = 'results',
    model: int = 0,
    tissuename: str = None,
    topgenenumber: int = 20,
    logfc_threshold: float = 0,
    add_info: str = None,
    username: str = None,
    password: str = None
) -> Union[str, List[str]]:
    """
    Annotate cell types using Llama 3 model via HuggingChat.
    
    Parameters:
    -----------
    input : Union[str, Dict[str, List[str]]]
        Either the path to a csv file or a dictionary of gene lists
    output_dir : str, default 'results'
        Directory to save the results
    model : int
        Index of the model to use.
    tissuename : str, optional
        Name of the tissue being analyzed
    topgenenumber : int, default 10
        Number of top differential genes to use if input is a csv
    logfc_threshold : float
        Log Fold-Change threshold for significance.
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
        # If input is a dictionary, nothing needs to be done
        print('Dictionnary input detected - make sure values are comma-separated gene lists')
        processed_input = input
        topgenenumber = list(test.values())[0].count(',') + 1
    elif isinstance(input, str):
        # If input is a string, filter significant genes and get top genes per cluster
        processed_input = seuratmarker_to_dict(input, 
                                               topgenenumber=topgenenumber, 
                                               logfc_threshold=logfc_threshold)
    else:
        raise ValueError("Input must be a string path to a csv or a dictionary of gene lists")
    
    # Construct prompt
    prompt = (
        f"Identify cell types of {tissuename} cell clusters using the following gene markers. "
        "Each row corresponds to one cluster. "
        "Only provide the cell type name. Do not show cluster numbers before the cell type name. "
        "Some could be a mixture of multiple cell types. "
        f"{add_info or ''}\n" +
        "\n".join([f"{cluster} : {genes}" for cluster, genes in processed_input.items()])
    )

    print(prompt)

    # Write the prompt to a file
    with open(f'{output_dir}/sbm_20_prompt.txt', 'w') as f:
        f.write(prompt)
    
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
        print('Note that other models are available and can be selected by changing the model number :\n')
        for i, mod in enumerate(available):
            print(f'{i} : {mod.name}')
        
        
        # Send prompt and get response   
        response = chatbot.chat(prompt).wait_until_done()
        annotations = response.split('\n')
        
        # Trim and validate annotations
        cell_type_annotations = [
            annotation.strip() 
            for annotation in annotations 
            if annotation.strip()     #and not annotation.startswith(tuple('0123456789'))
        ]
        
        # Ensure we have one annotation per cluster in the batch 
        # TODO
        
        
        print('\nAnnotation done !')
        print('ALWAYS check the results returned by this function in case of AI hallucination, before proceeding with downstream analysis.')
        
        model_name = available[model].name.split('/', 1)[-1]
        list_to_csv(cell_type_annotations, f'{output_dir}/{tissuename}_{model_name}_{topgenenumber}markergenes_annotation.csv')

        return cell_type_annotations
    
    except Exception as e:
        print(f"Error in HuggingChat interaction: {e}")
        return prompt


# Example usage
if __name__ == "__main__":
    
    # Test method
    file = '../data/raw_cp/cp_cluster_markers.csv'
    test = seuratmarker_to_dict(file,
                                topgenenumber=20)
    
    # Set HUGGINGCHAT_USERNAME and HUGGINGCHAT_PASSWORD environment variables before running
    # (This must be done independently of the script)
    annotations = huggingchatcelltype(
        input=test, 
        model=7,
        tissuename='mouse choroid plexus', 
        add_info='Be as precise as possible.'
    )

    print(annotations)



