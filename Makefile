# SHELL += -x

CXX = g++
CFLAG = "-std=c++11"

app: app_csr.cpp conf/config.hpp
	$(CXX) $(CFLAG) app_csr.cpp -o app -lasound -lpthread

.PHONY:clean
clean:
	rm app
