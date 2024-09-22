import json
import copy
import uuid
import ast
import os
import io
import requests
import re

# import pandas as pd
# import numpy as np
# import dask.array as da
import zarr
import duckdb
import tifffile
import s3fs
from anndata import AnnData

from termcolor import colored
from tenacity import retry, wait_random_exponential, stop_after_attempt
from dotenv import load_dotenv

from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain.chains import (
    create_retrieval_chain,
    create_history_aware_retriever,
    LLMChain,
    ConversationChain
)
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.pydantic_v1 import BaseModel, Field
from langchain_core.runnables import ConfigurableField
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain_community.chat_models.ollama import ChatOllama
from langchain_community.embeddings import GPT4AllEmbeddings
from langchain_community.document_loaders import WebBaseLoader, PyMuPDFLoader
from langchain_community.tools.tavily_search import TavilySearchResults
from langchain_text_splitters import RecursiveCharacterTextSplitter
# from langchain_chroma import Chroma
from langchain.prompts import PromptTemplate
from langchain.memory import ConversationBufferMemory
from ome_zarr.io import parse_url
import ome_types

from bs4 import BeautifulSoup as bs4

from vitessce import (
    VitessceConfig,
    Component as cm,
    CoordinationType as ct,
    CoordinationLevel as CL,
    OmeTiffWrapper,
    MultiImageWrapper,
    VitesscePlugin,
    AnnDataWrapper,
    get_initial_coordination_scope_prefix
)

# from deepmerge import always_merger as merge
# from esbuild_py import transform

import ipylangchat
import ipywidgets as widgets
from IPython.display import display, HTML

# from utils import *
from os.path import join
