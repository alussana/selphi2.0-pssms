FROM ubuntu
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y python3
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN apt-get install -y python3-pip
RUN pip3 install --break-system-packages h5py
RUN pip3 install --break-system-packages numpy
RUN pip3 install --break-system-packages pandas
RUN pip3 install --break-system-packages openpyxl
RUN pip3 install --break-system-packages requests
RUN pip3 install --break-system-packages scikit-learn
RUN pip3 install --break-system-packages scipy
RUN pip3 install --break-system-packages matplotlib
RUN pip3 install --break-system-packages seaborn
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/alussana/radialtree.git
RUN pip3 install --break-system-packages ./radialtree
