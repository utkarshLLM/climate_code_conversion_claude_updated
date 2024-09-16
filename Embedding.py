
import torch
########## Embedding ###############

# From Medium Article
# EMBEDDING_QPM = 100
# EMBEDDING_NUM_BATCH = 5
# embeddings = VertexAIEmbeddings(
#     requests_per_minute=EMBEDDING_QPM,
#     num_instances_per_batch=EMBEDDING_NUM_BATCH,
#     model_name = "textembedding-gecko@latest"
# )


########## Embedding ###############

# GPT-generated code.
from transformers import AutoModel, AutoTokenizer
from transformers import RobertaTokenizer, RobertaModel

def generate_embeddings(code_chunks):
    # Load a pre-trained model
    model_name1 = "microsoft/codebert-base" #Use the codeBERT
    model_name2 = "microsoft/graphcodebert-base" #Also captures the structural information, use when the whole embeddings and retrieval works.

    tokenizer = AutoTokenizer.from_pretrained(model_name1)
    model = AutoModel.from_pretrained(model_name1)

    # tokenizer = RobertaTokenizer.from_pretrained("microsoft/codebert-base")
    # model = RobertaModel.from_pretrained("microsoft/codebert-base")

    
    # Function to generate embeddings
    def get_embeddings(code_chunks, model, tokenizer, max_length=256):
        # Tokenize the code chunks
        ## these parameters can also be fine tuned: padding, truncation, return tensors
        inputs = tokenizer(code_chunks, padding=True, truncation=True, max_length=max_length, return_tensors="pt")
        
        # Generate embeddings
        with torch.no_grad():
            outputs = model(**inputs)
        
        # Use the [CLS] token embedding as the code chunk embedding
        embeddings = outputs.last_hidden_state[:, 0, :].numpy() ## CHECK THIS
        
        return embeddings


    # Example usage
    embeddings = get_embeddings(code_chunks, model, tokenizer)

    return embeddings




### Searching the weaviate:
# # Optional: Search Weaviate by embedding similarity
# query_embedding = all_embeddings[0]  # Use one of the embeddings for the search
# result = client.query.get(class_name).with_near_vector(query_embedding).with_limit(3).do()

# # Display search results
# print(json.dumps(result, indent=2))