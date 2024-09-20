FROM ubuntu
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y curl
RUN apt-get install -y python3
RUN rm /usr/lib/python*/EXTERNALLY-MANAGED
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN curl "https://bootstrap.pypa.io/get-pip.py" > get-pip.py
RUN python get-pip.py
RUN pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
RUN pip3 install transformers
RUN pip3 install jupyterlab
RUN pip3 install ipywidgets
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install scipy
COPY ./phosformer/ /phosformer/
WORKDIR phosformer