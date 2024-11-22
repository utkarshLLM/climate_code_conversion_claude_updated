from tree_sitter_script import generate_chunks
from Embedding import generate_embeddings
# from langchain_anthropic import ChatAnthropic
# from langchain_anthropic import AnthropicLLM
from langchain.vectorstores import FAISS
from langchain.schema import Document
import Vector_Database
from config import client
from config import model_name
from messages import system_translate
from fortran_code_ import fortran_code 

code_chunks = generate_chunks(fortran_code)
embeddings = generate_embeddings(code_chunks)

## VectorStore (LangChain)
# Create Index from embedded code chunks

# Create Document objects
documents = [Document(page_content=chunk) for chunk in code_chunks]

# Create a custom Embeddings class to use with FAISS
class CodeBERTEmbeddings:
    def embed_documents(self, texts):
        return generate_embeddings(texts)

    def embed_query(self, text):
        return generate_embeddings([text])[0]

# Create the FAISS index
embeddings_model = CodeBERTEmbeddings()
# Vector Store Index
db = FAISS.from_documents(documents, embeddings_model)


###### ^^ db Didn't Work ^^ ######

## Method 3.1: MetaData Retrieval
### Using llama_index with LangChain - hope its fine, lol
from llama_index.core.vector_stores import MetadataFilters, ExactMatchFilter

### How do you get the chunk_id ???

for code_chunk in code_chunks:

    specific_chunk_id = code_chunk.node_id # (????)

    filters = MetadataFilters(filters=[
        ExactMatchFilter(key="chunk_id", value=specific_chunk_id)
    ])

    query_engine = db.as_query_engine(filters=filters) ### Does this method work with db?
    response = query_engine.query("Retrieve all metadata for this chunk")

    print(response)

## Code llm
# code_llm = VertexAI(
#     model_name="code-bison@latest",
#     max_output_tokens=2048,
#     temperature=0.1,
#     verbose=False,
#     )

# Instantiate the Anthropic LLM with additional parameters
# through AnthropicLLM: Although this only supports Claude 2 legacy models
# code_llm = AnthropicLLM(
#     model=model_name,
#     max_tokens_to_sample=8092, ## max output
#     temperature=0, ## can tweak temperature for fine-tuning
#     top_p=1,
#     top_k=5,
#     streaming=False, 
# )

# Instantiate the Anthropic LLM with additional parameters
# through ChatAnthropic:
# code_llm = ChatAnthropic(
#     model=model_name,
#     max_tokens_to_sample=8092, ## max output
#     temperature=0, ## can tweak temperature for fine-tuning
#     top_p=1,
#     top_k=5,
#     streaming=False, 
# )

# How to get a response from the above^^
query = ""
# response = code_llm.invoke(query, )  ## max_tokens_to_sample is given here


## Query the vector database
### You have to make the query embedding too!
# query_embedding = get_embedding("! Another code snippet to query")
# result = vector_store.query([query_embedding], top_k=5) # top 5 
# print(result)

## create an LLM for generation


## see for integrating that into the larger framework for translation that we have.




