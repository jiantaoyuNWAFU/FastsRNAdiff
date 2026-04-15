from setuptools import setup,find_packages

setup(
        name="RegionRepCalc",
        version="0.1.0",
        py_modules=["RegionRepCalc"],
        entry_points={
            'console_scripts':[
                'RegionRepCalc=RegionRepCalc:main',
            ],
        },

        author="Xue S., Zhang X. , Gao Y. , etc.",
        author_email="jiantao.yu@nwafu.edu.cn",
        description="A tool to calculate total mapped reads and Rep-total reads for small RNA cluster",
    python_requires=">=3.7"
)
