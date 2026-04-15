from setuptools import setup, find_packages

setup(
    name="FastsRNAdiff",
    version="0.1.0",
    py_modules=["FastsRNAdiff", "Visualization"],
    install_requires=[
        "numpy >= 1.21.0",
        "pandas >= 1.5.0",
        "matplotlib >= 3.6.0",
        "seaborn >= 0.12.0",
        "tqdm >= 4.64.0",
        "scipy >= 1.9.0"
    ],
    entry_points={
        'console_scripts': [
            'FastsRNAdiff=FastsRNAdiff:main',
        ],
    },

    author="Xue S., Zhang X. , Gao Y. , etc.",
    author_email="jiantao.yu@nwafu.edu.cn",
    description="A highly efficient software for screening differentially expressed small RNA clusters",
    python_requires=">=3.7"
)
