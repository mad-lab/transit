From r-base:3.4.1
RUN apt-get update -y && apt-get install -y -f python2 python-dev python-pip
ADD src/ /src
ADD tests/ /tests
RUN pip install pytest 'numpy~=1.15' 'scipy~=1.2' 'matplotlib~=2.2' 'pillow~=5.0' 'statsmodels~=0.9' 'rpy2<2.9.0'
RUN R -e "install.packages('MASS')"
RUN R -e "install.packages('pscl')"

CMD [ "pytest", "./tests" ]

