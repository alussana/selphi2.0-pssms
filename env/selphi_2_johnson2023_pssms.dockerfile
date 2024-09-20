FROM ubuntu
RUN apt-get update
RUN apt-get upgrade -y 
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y python3
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN apt-get install -y python3-pip
RUN pip3 install h5py
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install openpyxl
RUN pip3 install requests
RUN pip3 install scikit-learn
RUN pip3 install scipy
RUN pip3 install matplotlib
RUN pip3 install seaborn
