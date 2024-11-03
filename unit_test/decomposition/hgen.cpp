#include <tfhe++.hpp>

template <class P>
void hGen() {
    constexpr std::array<typename P::T, P::l> hg = TFHEpp::hgen<P>();
    std::cout << std::dec << "l=" << P::l << std::endl;
    for(int i = 0; i < P::l; i++) {
        std::cout << std::hex  << " Ox"<< hg[i];
    }
    std::cout  << std::endl;

}

void hgen() {
    std::cout << "lvl1param"<< std::endl;
    hGen<TFHEpp::lvl1param>();
    std::cout << "lvl2param"<< std::endl;
    hGen<TFHEpp::lvl2param>();
    std::cout << "lvl3param"<< std::endl;
    hGen<TFHEpp::lvl3param>();
}


int main()
{
    hgen();
    return 0;
}
