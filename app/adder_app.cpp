#include "adder/adder.hpp"
#include <iostream>

int main(){
  int result = adder::add_one(1);
  std::cout << "1 + 1 = " << result << std::endl;
}
