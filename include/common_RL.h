#pragma once
#include <memory>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <sstream>
#include <fcntl.h>
#include <stdint.h>
#include "ReinforcementLearning.pb.h"
#include <armadillo>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <google/protobuf/message.h>

using google::protobuf::io::FileInputStream;
using google::protobuf::io::FileOutputStream;
using google::protobuf::io::ZeroCopyInputStream;
using google::protobuf::io::CodedInputStream;
using google::protobuf::io::ZeroCopyOutputStream;
using google::protobuf::io::CodedOutputStream;
using google::protobuf::Message;

namespace ReinforcementLearning{
inline bool ReadProtoFromTextFile(const char* filename, Message* proto) {
  int fd = open(filename, O_RDONLY);
//  CHECK_NE(fd, -1) << "File not found: " << filename;
  FileInputStream* input = new FileInputStream(fd);
  bool success = google::protobuf::TextFormat::Parse(input, proto);
  delete input;
//  close(fd);
  return success;
}


class RandomStream{
private:
    std::shared_ptr<std::mt19937> genPtr;
    std::shared_ptr<std::uniform_real_distribution<>> randomPtr_unitformReal;
    std::shared_ptr<std::uniform_int_distribution<>> randomPtr_unitformInt; 
public:     
    RandomStream(){
        
        std::random_device rd;
        genPtr = std::make_shared<std::mt19937>(rd());
        randomPtr_unitformReal = std::make_shared<std::uniform_real_distribution<>>(0.0, 1.0);
    }
    RandomStream(int low , int high){
        
        std::random_device rd;
        genPtr = std::make_shared<std::mt19937>(rd());
        
        randomPtr_unitformReal = std::make_shared<std::uniform_real_distribution<>>(0.0, 1.0);
        randomPtr_unitformInt = std::make_shared<std::uniform_int_distribution<>>(low, high);
    }
    double nextDou(){return (*randomPtr_unitformReal)(*genPtr);}
    int nextInt(){return (*randomPtr_unitformInt)(*genPtr);}
};



}
