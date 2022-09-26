from setuptools import setup, find_namespace_packages


setup(
    name='auto-cpo',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-cpo = auto_cpo.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description=' Automated analysis of carbapenemase-producing organism (CPO) sequence data',
    url='https://github.com/BCCDC-PHL/auto-cpo',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
