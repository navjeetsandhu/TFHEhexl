#include <tfhe++.hpp>

template <class P>
void totalOffset() {
    constexpr typename P::T offset = TFHEpp::totalOffset<P>();
    std::cout << std::hex << "Ox"<< std::hex << offset << std::endl;

}

void totaloffset() {
    std::cout << "lvl1param"<< std::endl;
    totalOffset<TFHEpp::lvl1param>();
    std::cout << "lvl2param"<< std::endl;
    totalOffset<TFHEpp::lvl2param>();
    std::cout << "lvl3param"<< std::endl;
    totalOffset<TFHEpp::lvl3param>();
}


int main()
{
    totaloffset();
    return 0;
}
