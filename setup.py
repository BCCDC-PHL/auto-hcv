from setuptools import setup, find_namespace_packages


setup(
    name='auto-hcv',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-hcv = auto_hcv.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description=' Automated analysis of HCV',
    url='https://github.com/BCCDC-PHL/auto-hcv',
    author='Sherrie Wang, Dan Fornika',
    author_email='sherrie.wang@bccdc.ca, dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
