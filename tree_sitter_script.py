import tree_sitter_fortran as tsfortran
from tree_sitter import Language, Parser
from tree_sitter_languages import get_language, get_parser
from llama_index.core.node_parser import CodeSplitter


def generate_chunks(fortran_code: str):
    FORTRAN_LANGUAGE = Language(tsfortran.language())
    parser = Parser(FORTRAN_LANGUAGE)

    language = "fortran"

    splitter_fortran = CodeSplitter(
        language = language,
        parser = parser,
        chunk_lines = 50,  # lines per chunk
        chunk_lines_overlap = 15,  # lines overlap between chunks
        max_chars = 1500  # max chars per chunk
    )

    split_code = splitter_fortran.split_text(fortran_code)

    return split_code

    ### This for checking if the chunks are working or not
    # with open("/Users/utkarshprajapati/Desktop/Columbia Course Files/RA - Climate Modelling/climate_code_conversion/_codegen_claude/micro_mg1_0.F90", "r") as file:
    #     fortran_code = file.read()