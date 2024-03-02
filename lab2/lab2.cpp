#include <stdint.h>
#include <vector>
#include <iostream>
#include <bitset>

// регистр на 16 бит
void LFSR16(std::vector<uint8_t>& coefficients)
{
	uint16_t startOutput = 0xACE1u;
	uint16_t currentOutput = startOutput;
	uint16_t tmp;
	unsigned t = 0;
	do
	{
		tmp = 0;
		for (auto i = coefficients.begin(); i != coefficients.end(); ++i)
		{
			tmp = tmp ^ (currentOutput >> (15 - *i));
		}
		currentOutput = (currentOutput >> 1) | (tmp << (15));
		std::bitset<16> bits(currentOutput);
		std::cout << t << ":" << currentOutput << std::endl;
		t++;
	} while (currentOutput != startOutput);
	std::cout << "Perior = " << t << std::endl;
}

// регистр для n <= 32 вроде работает, но по размеру числа все равно будут 4 байта
void LFSR(std::vector<uint8_t>& coefficients, uint32_t startValue)
{
	uint8_t n = coefficients[0];
	uint32_t startOutput = startValue;
	uint32_t currentOutput = startOutput;
	uint32_t tmp;
	uint32_t mask = ((1 << (n)) - 1);
	std::bitset<32> mask_bits(mask);
	unsigned t = 0;
	do
	{
		tmp = 0;
		for (auto i = coefficients.begin(); i != coefficients.end(); ++i)
		{
			tmp = tmp ^ (currentOutput >> (n-*i));
		}
		currentOutput = (((currentOutput >> 1) | (tmp << (n-1))) & mask);
		std::bitset<32> bits(currentOutput);
		std::cout << t << ":" << currentOutput << std::endl;
		t++;
	} 
	while (currentOutput != startOutput);
	std::cout << "Perior = " << t << std::endl;
}

int main() {
	//std::vector<uint8_t> coeff = { 16, 14, 13, 11};
	//LFSR16(coeff);
	std::vector<uint8_t> coeff = { 2, 1};
	LFSR(coeff, 1);
	return 0;
}