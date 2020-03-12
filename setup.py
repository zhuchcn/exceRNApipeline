from setuptools import setup, find_packages, Extension
import sys

extra_link_args = ['-lz']

if sys.platform == 'darwin':
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = "10.9"
    extra_link_args += ['-bundle', '-flat_namespace', '-undefined', 'suppress']

CONSOLE_SCRIPTS = [
    'pipeline=exceRNApipeline.pipeline.__main__:main',
    'task_process.py=exceRNApipeline.tasks.task_preprocess:main',
    'task_preprocess_stats.py=exceRNApipeline.tasks.task_preprocess_stats:main',
    'task_star_index.py=exceRNApipeline.tasks.task_star_index:main',
    'task_star_align.py=exceRNApipeline.tasks.task_star_align:main',
    'task_anno_to_fasta.py=exceRNApipeline.tasks.task_anno_to_fasta:main',
    'task_solve_silva_taxa.py=exceRNApipeline.tasks.task_solve_silva_taxa:main',
    'task_extract_fasta.py=exceRNApipeline.tasks.task_extract_fasta:main',
    'task_summarize_counts.py=exceRNApipeline.tasks.task_summarize_counts:main',
    'task_endogenous_count_barplot.py=' + \
            'exceRNApipeline.tasks.task_endogenous_count_barplot:main',
    'task_silva_extract_unmapped.py=exceRNApipeline.tasks.task_silva_extract_unmapped:main',
    'task_split_fasta.py=exceRNApipeline.tasks.task_split_fasta:main',
    'task_bacteria_download_genome.py=' + \
            'exceRNApipeline.tasks.task_bacteria_download_genome:main',
    'task_star_align_bacteria.py=exceRNApipeline.tasks.task_star_align_bacteria:main',
    'task_count_exogenous_taxa.py=exceRNApipeline.tasks.task_count_exogenous_taxa:main'
]

MODULE_TAXA = Extension(
    name="_exceRNApipeline_taxaCounter",
    sources=[
        'exceRNApipeline/includes/taxonomy/taxonomy.i',
        'exceRNApipeline/includes/taxonomy/NCBITaxonomy.cpp',
        'exceRNApipeline/includes/taxonomy/TaxaCounter.cpp',
        'exceRNApipeline/includes/gzstream/gzstream.C'
    ],
    swig_opts=['-c++'],
    extra_compile_args=['-IexceRNApipeline/includes/gzstream', '-lz', '-g', '-fPIC'],
    extra_link_args=extra_link_args
)

setup(
    name = 'exceRNApipeline',
    version = '0.1.0',
    author = 'Chenghao Zhu',
    author_email = 'zhuchcn@gmail.com',
    install_requires = [
        "snakemake",
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        'pathos',
        'python-magic',
        'biopython'
    ],
    ext_modules=[MODULE_TAXA],
    packages=find_packages(exclude=["exceRNApipeline/tests*"]),
    package_data={"exceRNApipeline.pipeline": ["smk/*"]},
    entry_points = {
        'console_scripts': CONSOLE_SCRIPTS
    }
)