#include "c_assert.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <string>
#include <sstream>


using namespace std;
using namespace TFHEpp;

constexpr int batch = otherparam::batch;
vector<uint8_t> pa(batch);
vector<uint8_t> pb(batch);
vector<uint8_t> pres(batch);

TLWEn<lvl1param, batch> ca;
TLWEn<lvl1param, batch> cb;
TLWEn<lvl1param, batch> cres;


int main()
{

    cout << "batch: " << batch << endl;
    int j;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);

    for (j = 0; j < batch; j++) pa[j] = binary(engine) > 0;
    for (j = 0; j < batch; j++) pb[j] = binary(engine) > 0;

    vector<TLWE<TFHEpp::lvl1param>> caa(batch);
    vector<TLWE<TFHEpp::lvl1param>> cbb(batch);
    vector<TLWE<TFHEpp::lvl1param>> ccres(batch);

    caa = bootsSymEncrypt(pa, *sk);
    cbb = bootsSymEncrypt(pb, *sk);

    for (j = 0; j < batch; j++) {
        ca[j] = caa[j];
        cb[j] = cbb[j];
    }

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    HomNANDbatch(cres, ca, cb, ek);

    for (j = 0; j < batch; j++) {
        ccres[j] = cres[j];
    }

    end = chrono::system_clock::now();

    pres = bootsSymDecrypt(ccres, *sk);
    for (int i = 0; i < batch; i++) {
        c_assert(pres[i] == !(pa[i] & pb[i]));
    }
   cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / batch << "ms" << endl;
}