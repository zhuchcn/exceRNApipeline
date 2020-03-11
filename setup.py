from setuptools import setup, find_packages, Extension
import sys

extra_link_args = ['-lz']

if sys.platform == 'darwin':
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = "10.9"
    extra_link_args += ['-bundle', '-flat_namespace', '-undefined', 'suppress']

CONSOLE_SCRIPTS = [
    'task_process.py=src.tasks.task_preprocess:main',
    'task_preprocess_stats.py=src.tasks.task_preprocess_stats:main',
    'task_star_index.py=src.tasks.task_star_index:main',
    'task_star_align.py=src.tasks.task_star_align:main',
    'task_anno_to_fasta.py=src.tasks.task_anno_to_fasta:main'
    'task_solve_silva_taxa.py=src.tasks.task_solve_silva_taxa:main',
    'task_extract_fasta.py=src.tasks.task_extract_fasta:main',
    'task_summarize_counts.py=src.tasks.task_summarize_counts.py'
    'task_endogenous_count_barplot.py=' + \
            'src.tasks.task_endogenous_count_barplot.py',
    'task_silva_extract_unmapped.py=src.tasks.task_silva_extract_unmapped.py',
    'task_split_fasta.py=src.tasks.task_split_fasta.py',
    'task_bacteria_download_genome.py=' + \
            'src.tasks.task_bacteria_download_genome.py'
    'task_star_align_bacteria.py=src.tasks.task_star_align_bacteria.py'
]

module_taxa = Extension(
    name="_exRNApipeline_taxaCounter",
    sources=[
        'src/includes/taxonomy/taxonomy.i',
        'src/includes/taxonomy/NCBITaxonomy.cpp',
        'src/includes/taxonomy/TaxaCounter.cpp',
        'src/includes/gzstream/gzstream.C'
    ],
    swig_opts=['-c++'],
    extra_compile_args=['-Isrc/includes/gzstream', '-lz', '-g', '-fPIC'],
    extra_link_args=extra_link_args
)

setup(
    name = 'exRNApipeline',
    version = '0.1.0',
    author = 'Chenghao Zhu',
    author_email = 'zhuchcn@gmail.com',
    install_requires = [
        "snakemake",
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn"
    ],
    ext_modules=[module_taxa],
    packages = find_packages(exclude=["src.tests*"]),
    entry_points = {
        'console_scripts': CONSOLE_SCRIPTS
    }
)