from langchain.vectorstores import Weaviate
# import weaviate
# from weaviate.embedded import EmbeddedOptions


# def create_weaviate_client():
#     # Create a Weaviate client
#     client = weaviate.Client(
#         embedded_options=EmbeddedOptions()
#     )
#     return client

# def vector_store(code_chunks, embeddings):
#     # Define the class for our code snippets
#     class_name = "CodeSnippet"
    
#     client = create_weaviate_client()
#     # Check if the class already exists, if not create it
#     if not client.schema.exists(class_name):
#         class_obj = {
#             "class": class_name,
#             "vectorizer": "none",  # We'll provide our own vectors
#             "properties": [
#                 {
#                     "name": "content",
#                     "dataType": ["text"],
#                 }
#             ]
#         }
#         client.schema.create_class(class_obj)
    
#     # Create a Weaviate vector store
#     vector_store = Weaviate(
#         client=client,
#         index_name=class_name,
#         text_key="content",
#         embedding=embeddings,
#         by_text=False
#     )
    
#     # Add data to the vector store
#     vector_store.add_texts(code_chunks)
    
#     return vector_store

