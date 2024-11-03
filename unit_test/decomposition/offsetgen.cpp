#include <tfhe++.hpp>


template <class P>
void offsetGen() {
    constexpr typename P::T offset = TFHEpp::offsetgen<P>();
    std::cout << std::hex << "Ox"<< std::hex << offset << std::endl;
}

int main()
{
    std::cout << "lvl1param"<< std::endl;
    offsetGen<TFHEpp::lvl1param>();
    std::cout << "lvl2param"<< std::endl;
    offsetGen<TFHEpp::lvl2param>();
    std::cout << "lvl3param"<< std::endl;
    offsetGen<TFHEpp::lvl3param>();

    return 0;
}
