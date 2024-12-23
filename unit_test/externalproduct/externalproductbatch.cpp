#include "c_assert.hpp"
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    constexpr int batch = 27;
    cout << "test p=1: lvl1 batch" << endl;

    {
        lweKey key;
        BooleanArrayn<lvl1param::n, batch> p;

        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++)
                p[j][i] = (binary(engine) > 0);

        Polynomialn<lvl1param, batch> pmu;
        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++)
                pmu[j][i] = p[j][i] ? lvl1param::mu : -lvl1param::mu;

        TRLWEn<lvl1param, batch> c = trlweSymEncryptbatch<lvl1param, batch>(pmu, key.lvl1);

        Polynomialn<TFHEpp::lvl1param, batch> plainpoly = {
            static_cast<typename lvl1param::T>(0)};

        for (int j = 0; j < batch; j++)
            plainpoly[j][0] = 1;

        //for (int j = 0; j < batch; j++)
        //    for (int i = 0; i < lvl1param::n; i++)
         //       cout << j << " " << i << " " << plainpoly[j][i] <<  endl;

        TRGSWFFTn<lvl1param, batch> trgswfft =
            trgswfftSymEncryptbatch<lvl1param, batch>(plainpoly, key.lvl1);
        trgswfftExternalProductbatch<lvl1param, batch>(c, c, trgswfft);


        BooleanArrayn<lvl1param::n, batch> p2 = trlweSymDecryptbatch<lvl1param, batch>(c, key.lvl1);
        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++) {
                //cout << j << " " << i << " " << p[j][i] << "  " << p2[j][i] << endl;
                c_assert(p[j][i] == p2[j][i]);
            }
    }
    cout << "Passed" << endl;
    cout << "test p=-1: lvl1 batch" << endl;

    {
        lweKey key;
        BooleanArrayn<lvl1param::n, batch> p;

        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++)
                p[j][i] = (binary(engine) > 0);

        Polynomialn<lvl1param, batch> pmu;
        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++)
                pmu[j][i] = p[j][i] ? lvl1param::mu : -lvl1param::mu;

        TRLWEn<lvl1param, batch> c = trlweSymEncryptbatch<lvl1param, batch>(pmu, key.lvl1);

        Polynomialn<TFHEpp::lvl1param, batch> plainpoly = {
            static_cast<typename lvl1param::T>(0)};

        for (int j = 0; j < batch; j++)
            plainpoly[j][0] = -1;

        //for (int j = 0; j < batch; j++)
        //    for (int i = 0; i < lvl1param::n; i++)
        //       cout << j << " " << i << " " << plainpoly[j][i] <<  endl;

        TRGSWFFTn<lvl1param, batch> trgswfft =
            trgswfftSymEncryptbatch<lvl1param, batch>(plainpoly, key.lvl1);
        trgswfftExternalProductbatch<lvl1param, batch>(c, c, trgswfft);


        BooleanArrayn<lvl1param::n, batch> p2 = trlweSymDecryptbatch<lvl1param, batch>(c, key.lvl1);
        for (int j = 0; j < batch; j++)
            for (int i = 0; i < lvl1param::n; i++) {
                //cout << j << " " << i << " " << p[j][i] << "  " << p2[j][i] << endl;
                c_assert(p[j][i] == !p2[j][i]);
            }
    }
    cout << "Passed" << endl;


}