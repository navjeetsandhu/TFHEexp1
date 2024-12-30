#include "c_assert.hpp"
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <chrono>

using namespace std;
using namespace TFHEpp;

int main()
{

    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    cout << "lvl1 test p=1" << endl;
    chrono::system_clock::time_point start, end;
	double elapsed;

	start = chrono::system_clock::now();
    for (int testi = 0; testi< 10; testi++) {
        lweKey key;

        array<bool, lvl1param::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        Polynomial<lvl1param> pmu;
        for (int i = 0; i < lvl1param::n; i++)
            pmu[i] = p[i] ? lvl1param::mu : -lvl1param::mu;
        TRLWE<lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);

        const Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename lvl1param::T>(1)};

        //for (int i = 0; i < lvl1param::n; i++)
        //    cout << i << " " << plainpoly[i] <<  endl;


        TRGSWFFT<lvl1param> trgswfft =
            trgswfftSymEncrypt<lvl1param>(plainpoly, key.lvl1);
        trgswfftExternalProduct<lvl1param>(c, c, trgswfft);
        array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) {
            //cout << p[i] << "  " << p2[i] << endl;
            c_assert(p[i] == p2[i]);
        }
    }
	end = chrono::system_clock::now();

    elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / 10 << "ms" << endl;

    cout << "Passed" << endl;

    cout << "lvl1 test p=-1" << endl;

	start = chrono::system_clock::now();
    {
        lweKey key;

        array<bool, lvl1param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<typename TFHEpp::lvl1param::T, lvl1param::n> pmu;
        for (int i = 0; i < lvl1param::n; i++)
            pmu[i] = p[i] ? lvl1param::mu : -lvl1param::mu;
        TRLWE<lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);

        const Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename lvl1param::T>(-1)};

        TRGSWFFT<lvl1param> trgswfft =
            trgswfftSymEncrypt<lvl1param>(plainpoly, key.lvl1);
        trgswfftExternalProduct<lvl1param>(c, c, trgswfft);
        array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) {
            //cout << p[i] << "  " << p2[i] << endl;
            c_assert(p[i] == !p2[i]);
        }
    }

	end = chrono::system_clock::now();

    elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed  << "ms" << endl;

    cout << "Passed" << endl;

}