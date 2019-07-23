From r-base:3.6.1
RUN apt-get update -y && apt-get install -y -f python3 python-dev python3-pip
ADD src/ /src
ADD tests/ /tests
RUN pip3 install pytest 'numpy~=1.16' 'scipy~=1.2' 'matplotlib~=3.0' 'pillow~=6.0' 'statsmodels~=0.9' 'rpy2'
RUN R -e "install.packages('MASS')"
RUN R -e "install.packages('pscl')"

CMD [ "pytest", "./tests" ]

