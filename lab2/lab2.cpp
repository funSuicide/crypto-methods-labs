#include <stdint.h>
#include <vector>
#include <iostream>
#include <bitset>
#include <fstream>

uint32_t moduleTransform(uint32_t value, uint8_t p)
{
	return (value % p + p) % p;
}

// регистр для n <= 32 вроде работает, но по размеру числа все равно будут 4 байта
std::vector<uint8_t> LFSR(std::vector<uint8_t>& coefficients, uint32_t startValue)
{
	std::vector<uint8_t> result;

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

		uint8_t currentBit = (currentOutput >> (n-1)) & 1u;
		result.push_back(currentBit);
		t++;
	} 
	while (currentOutput != startOutput);
	//std::cout << "Period = " << t << std::endl;

	return result;
}

std::vector<uint8_t> generator(uint32_t size)
{
	std::vector<uint8_t> coeff = { 20, 17 };
	uint8_t startValue = 0xA;

	std::vector<uint8_t> result(size);
	for (uint32_t currentSize = 0; currentSize < size; )
	{
		std::vector<uint8_t> x = LFSR(coeff, 0xA2);
		std::vector<uint8_t> a = LFSR(coeff, 0xB7);

		for (size_t i = 0; i < a.size(); ++i)
		{
			if (a[i] == 1)
			{
				result[currentSize] = x[i];
				currentSize++;
				std::cout << currentSize << std::endl;
				if (currentSize >= size)
				{
					break;
				}
			}
		}
	}
	return result;
}

int main() {
	
	uint32_t size = 1000000;

	std::vector<uint8_t> bits = generator(size);

	std::ofstream file("test.txt", std::ios_base::out);
	if (file.is_open()) { 
		for (size_t i = 0; i < bits.size(); ++i)
		{
			if (bits[i] == 0)
			{
				file << '0';
			}
			else
			{
				file << '1';
			}
		}
		file.close();
	}


	
	return 0;
}