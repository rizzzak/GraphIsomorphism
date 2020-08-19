//Программа для поиска изоморфного подграфа в графе эвристическими методами.
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include <ctime>
#include <set>
#include <algorithm>
#include <string>
#include <iterator>
#include <fstream>
#include <thread>
using namespace std;

int iterCnt = 200;
double epsilon = 0.0001;


bool NextSet(vector<int> &a, int n, int n2)
{
	for (int i = n - 1; i >= 0; --i)
		if (a[i] < n2 - n + i)
		{
			++a[i];
			for (int j = i + 1; j < n; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}
void printVectorInt(vector<int> vec, string name)
{
	cout << "\n" << name << ", размер = " << vec.size() << ": ";
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << " ";
}
void printMatrInt(vector<vector<int>> matr, string name)
{
	cout << "\n" << name << ", размер = " << matr.size() << ": \n";
	for (int i = 0; i < matr.size(); i++)
	{
		for (int j = 0; j < matr[i].size(); j++)
		{
			cout << matr[i][j] << " ";
		}
		cout << "\n";
	}
}
void printVectorOfMatrInt(vector<vector<vector<int>>> vec, string name)
{
	for (int k = 0; k < vec.size(); k++)
	{
		cout << "\nMatrix " << k << " - " << name << ", size = " << vec[k].size() << ": \n";
		for (int i = 0; i < vec[k].size(); i++)
		{
			for (int j = 0; j < vec[k][i].size(); j++)
			{
				cout << vec[k][i][j] << " ";
			}
			cout << "\n";
		}
		cout << "\n";
	}
}
int compareTables(vector<vector<int>> table, vector<vector<int>> table2, int n)
{
	int q = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (table[i][j] != table2[i][j])
				q++;
		}
	}
	return q;
}
vector<vector<int>> buildGraphFromP(vector<int> p, vector<vector<int>> source)
{
	vector<vector<int>> graph;
	vector<int> graphVertix;
	for (int i = 0; i < p.size(); i++)
	{
		graphVertix.push_back(0);
	}
	for (int i = 0; i < p.size(); i++)
	{
		for (int j = 0; j < p.size(); j++)
		{
			graphVertix[j] = source[ p[i] ] [ p[j] ];
		}
		graph.push_back(graphVertix);
	}
	return graph;
}
int swapRandomIteration(vector<vector<int>> table, vector<vector<int>> table2, vector<int>& perestanovka)
{
	int q, randV, randV2, n = table2.size(), n2 = table.size();
	vector<int> shortPerestanovka;
	while (true)
	{
		randV = rand() % n2;
		randV2 = rand() % n2;
		if (randV != randV2)
			break;
	}
	swap(perestanovka[randV], perestanovka[randV2]);
	for (int i = 0; i < n; i++)
		shortPerestanovka.push_back(perestanovka[i]);
	table = buildGraphFromP(shortPerestanovka, table);
	q = compareTables(table, table2, n);
	return q;	
}
int swapRandomWeightedIteration(vector<vector<int>> table, vector<vector<int>> table2, vector<int>& perestanovka, double d, int startQ)
{
	int quality,minF, newF, index, n = table2.size(), n2 = table.size();
	vector<int> currentSwap, shortPerestanovka, qualities, oldPerestanovka = perestanovka;
	vector<vector<int>> allSwaps, newTable;
	for (int i = 0; i < n2 - 1; i++)
	{
		for (int j = i + 1; j < n2; j++)
		{
			swap(perestanovka[i], perestanovka[j]);
			for (int k = 0; k < n; k++)
				shortPerestanovka.push_back(perestanovka[k]);
			newTable = buildGraphFromP(shortPerestanovka, table);
			quality = compareTables(newTable, table2, n);
			qualities.push_back(quality);
			perestanovka = oldPerestanovka;
			currentSwap.push_back(i);
			currentSwap.push_back(j);
			allSwaps.push_back(currentSwap);
			currentSwap.clear();
			shortPerestanovka.clear();
		}
	}
	double rk = (rand() % 10000) / 10000.;
	minF = (qualities[0]-startQ) * (1 + 2 * d * (rk - 0.5));
	index = 0;
	for (int i = 1; i < qualities.size(); i++)
	{
		rk = (rand() % 10000) / 10000.;
		newF = (qualities[i] - startQ) * (1 + 2 * d * (rk - 0.5));
		if (newF < minF)
		{
			minF = newF; 
			index = i;
		}
	}
	currentSwap = allSwaps[index];
	swap(perestanovka[allSwaps[index][0]], perestanovka[allSwaps[index][1]]);
	return qualities[index];
}
int swapGreedIteration(vector<vector<int>> table, vector<vector<int>> table2, vector<int>& perestanovka,int startQ)
{
	int n = table2.size();
	int n2 = table.size();
	int index;
	vector<int> qualities;
	int quality, minF, newF;
	int bestQuality, bestDeltaQuality(n * n), newDeltaQuality;
	vector<int> bestPerestanovka = perestanovka;
	vector<int> oldPerestanovka = perestanovka;
	vector<int> shortPerestanovka;
	vector<int> currentSwap;
	vector<vector<int>> allSwaps;
	vector<vector<int>> newTable;
	bestQuality = startQ;
	for (int i = 0; i < n2 - 1; i++)
	{
		for (int j = i + 1; j < n2; j++)
		{
			if (i == j)
				continue;
			swap(perestanovka[i], perestanovka[j]);
			for (int k = 0; k < n; k++)
				shortPerestanovka.push_back(perestanovka[k]);
			newTable = buildGraphFromP(shortPerestanovka, table);
			quality = compareTables(newTable, table2, n);
			newDeltaQuality = quality - startQ;
			if (newDeltaQuality < bestDeltaQuality)
			{
				bestDeltaQuality = newDeltaQuality;
				bestPerestanovka = perestanovka;
				bestQuality = quality;
			}
			perestanovka = oldPerestanovka;
			shortPerestanovka.clear();
		}
	}
	perestanovka = bestPerestanovka;
	return bestQuality;
}
int swapFullIteration(vector<vector<int>> table, vector<vector<int>> table2, vector<int>& perestanovka)
{
	int q, n = table2.size(), n2 = table.size();
	vector<int> shortPerestanovka;
	while (NextSet(perestanovka, n, n2))
	{
		for (int i = 0; i < n; i++)
			shortPerestanovka.push_back(perestanovka[i]);
		table = buildGraphFromP(shortPerestanovka, table);
		q = compareTables(table, table2, n);
		return q;
	}
	return (n * n - n);
}
vector<vector<int>> genGraphWithDensity(int n, double D)
{
	vector<vector<int>> newGraph;
	for (int i = 0; i < n; i++)
	{
		newGraph.push_back({});
		for (int j = 0; j < n; j++)
		{
			newGraph[i].push_back(0);
		}
	}
	double stepDPerEdge = 2. / (n * n - n);
	int edgesCnt;
	for (int k = 0; k < ((n * n - n) / 2.) + 1.0; k++)
	{
		if (k != 0)
		{
			if ((D > k * stepDPerEdge - epsilon - stepDPerEdge / 2) && (D <= k * stepDPerEdge + epsilon + stepDPerEdge / 2))
			{
				edgesCnt = k;
				break;
			}
		}
		else
		{
			if (D <= stepDPerEdge / 2)
			{
				edgesCnt = k;
				break;
			}
		}
	}
	int randomVertix1;
	int randomVertix2;
	for (int i = 0; i < edgesCnt; i++)
	{
		while (true)
		{
			randomVertix1 = rand() % n;
			randomVertix2 = rand() % n;
			if (randomVertix1 == randomVertix2)
				continue;
			if (newGraph[randomVertix1][randomVertix2] == 0)
			{
				newGraph[randomVertix1][randomVertix2] = 1;
				newGraph[randomVertix2][randomVertix1] = 1;
				break;
			}
		}
	}
	return newGraph;
}
long double fact(int N)
{
	if (N < 0) // если пользователь ввел отрицательное число
		return 0; // возвращаем ноль
	if (N == 0) // если пользователь ввел ноль,
		return 1; // возвращаем факториал от нуля - не удивляетесь, но это 1 =)
	else // Во всех остальных случаях
		return N * fact(N - 1); // делаем рекурсию.
}
int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Добро пожаловать в программу для поиска подграфа в графе G изоморфного графу G'.\n";
	cout << "Версия: 1.0.\n";
	srand(time(NULL));
	vector<vector<int>> table = { { 0,1,0,0,1 },{ 1,0,1,1,1 },{ 0,1,0,0,1 },{ 0,1,0,0,0 },{ 1,1,1,0,0 } };
	vector<vector<int>> table2 = { { 0,0,1,0 },{ 0,0,1,1 },{ 1,1,0,1 },{ 0,1,1,0 } };
	vector<int> perestanovka;
	//for sample
	double averageQRandom(0.), averageQGreed(0.), averageQRW(0.), averageQFull(0.);
	vector<double> QRandom, QGreed, QRW, QFull;
	double totalProbabilityRandom(0.), totalProbabilityGreed(0.), totalProbabilityRW(0.), totalProbabilityFull(0.);
	double dispersionRandom(0.), dispersionGreed(0.), dispersionRW(0.), dispersionFull(0.);
	double sigmaRandom, sigmaGreed, sigmaRW, sigmaFull;
	double deltaQRandom, deltaQGreed, deltaQRW, deltaQFull;
	//for subgraphs
	double localAverageQRandom(0.), localAverageQGreed(0.), localAverageQRW(0.), localAverageQFull(0.);
	int localQRandom, localQGreed, localQRW, localQFull;
	double localProbabilityRandom(0.), localProbabilityGreed(0.), localProbabilityRW(0.), localProbabilityFull(0.);
	int localMinQRandom, localMinQGreed, localMinQRW, localMinQFull;
	int skippedSubGraphs(0);
	bool mute = true;
	double density;
	int n2, n;
	double d = 0.1;
	cout << "Данная программа позволяет решать 2 задачи:\n";
	cout << " 1 - Оценка эффективности работы эвристических алгоритмов для поставленной задачи на случайной выборке графов\n";
	cout << " 2 - Поиск изоморфного подграфа в заданном графе с помощью эвристических алгоритмов\n";
	int type = 0;
	cout << "\nПожалуйста, введите тип работы программы(1 или 2): ";
	cin >> type;
	while ((type != 1) && (type != 2))
	{
		cout << "Ошибка 1: Введен неверный режим работы. Повторите ввод.";
		cin >> type;
	}
	if (type == 1)
	{
		ofstream outRWMeta;
		ofstream outRandom, outRW, outGreed, outFull;
		outFull.open("Full.csv");
		outRW.open("RW.csv");
		outRandom.open("Random.csv");
		outGreed.open("Greed.csv");
		outRWMeta.open("RWMeta.csv");
		//////////////////////////////////////////////////////////
		//			Ф О Р М И Р У Е М	В Ы Б О Р К У
		vector<vector<vector<int>>> sampleTable, sampleTable2;
		int sizeOfSample = 2;
		n = 10;
		n2 = 20;
		cout << "\nВыбран первый тип работы - работа с выборкой.\n";
		cout << "Введите размерности графов: большого - n2(от 3 до 99) и малого - n(от 3 до 99).\n";
		cout << "n = ";
		cin >> n;
		while (!((n > 2) && (n < 100)))
		{
			cout << "Ошибка 7. Превышение диапазона данных переменной n. Повторите ввод.\n";
			cout << "n = ";
			cin >> n;
		}
		
		cout << "n2 = ";
		cin >> n2;
		while (!((n2 > 2) && (n2 < 100)))
		{
			cout << "Ошибка 7. Превышение диапазона данных переменной n2. Повторите ввод.\n";
			cout << "n2 = ";
			cin >> n2;
		}
		while (n2 <= n)
		{
			cout << "Ошибка 8. Размерность переменной n больше или равна размерности n2. Повторите ввод.\n";
			cout << "n = ";
			cin >> n;
			while (!((n > 2) && (n < 100)))
			{
				cout << "Ошибка 7. Превышение диапазона данных переменной n. Повторите ввод.\n";
				cout << "n = ";
				cin >> n;
			}
		
			cout << "n2 = ";
			cin >> n2;
			while (!((n2 > 2) && (n2 < 100)))
			{
				cout << "Ошибка 7. Превышение диапазона данных переменной n2. Повторите ввод.\n";
				cout << "n2 = ";
				cin >> n2;
			}
		}
		cout << "\nЧисло различных подграфов при заданных n2 и n: " << (fact(n2) / fact(n) / fact(n2 - n)) << "\n";
		cout << "Введите размер выборки - количество пар графов, на которых будет проводиться тест(3-999): ";
		cin >> sizeOfSample;
		while (!((sizeOfSample > 2) && (sizeOfSample < 1000)))
		{
			cout << "Ошибка 7. Превышение диапазона данных переменной sizeOfSample. Повторите ввод.\n";
			cout << "sizeOfSample = ";
			cin >> sizeOfSample;
		}
		cout << "Введите число итераций каждого из эвристических алгоритмов(3-999): ";
		cin >> iterCnt;
		while (!((iterCnt > 2) && (iterCnt < 1000)))
		{
			cout << "Ошибка 7. Превышение диапазона данных переменной iterCnt. Повторите ввод.\n";
			cout << "iterCnt = ";
			cin >> iterCnt;
		}
		cout << "\n";

		outRWMeta << "iterCnt = " << to_string(iterCnt) << "\n";
		outRWMeta << "sizeOfSample = " << to_string(sizeOfSample) << "\n";

		cout << "\nn = " << to_string(n) << ", n2 = " << to_string(n2);
		for (density = 0.2; density < 1.0; density += 0.2)
		{
			cout << "\n\nАнализ плотности графа " << density;
			cout << "\nВыполняется метаоптимизация метода случайного взвешенного перебора:\n";
			outFull << "Delta = " << to_string(n2 - n) << " D = " << to_string(density) << ";";
			outRW << "Delta = " << to_string(n2 - n) << " D = " << to_string(density) << ";";
			outRandom << "Delta = " << to_string(n2 - n) << " D = " << to_string(density) << ";";
			outGreed << "Delta = " << to_string(n2 - n) << " D = " << to_string(density) << ";";
			for (int i = 0; i < sizeOfSample; i++)
			{
				table = genGraphWithDensity(n2, density);
				table2 = genGraphWithDensity(n, density);
				sampleTable.push_back(table);
				sampleTable2.push_back(table2);
			}
			vector<int> pRW(perestanovka);
			d = 0.15;
			vector<double> totalMinQRWMO;
			outRWMeta << "n = " << to_string(n) << " n2 = " << to_string(n2) << " D = " << to_string(density) << ",";
			for (d = 0; d < 1.01; d += 0.05)
			{
				localQRW = n * n;
				localMinQRW = n * n;
				for (int j = 0; j < n2; j++)
				{
					perestanovka.push_back(j);
				}
				pRW = perestanovka;
				int startQ = compareTables(sampleTable[0], sampleTable2[0], n);
				for (int iteration = 0; iteration < iterCnt; iteration++)
				{
					localQRW = swapRandomWeightedIteration(sampleTable[0], sampleTable2[0], pRW, d, startQ);
					if (localMinQRW > localQRW)
						localMinQRW = localQRW;
				}
				totalMinQRWMO.push_back(localMinQRW);
				perestanovka.clear();
				cout << "x";
				outRWMeta << to_string(localMinQRW) << ",";
			}
			int minElementIndex = 0;
			localMinQRW = n * n;
			for (int i = 0; i < totalMinQRWMO.size(); i++)
			{
				if (localMinQRW > totalMinQRWMO[i])
				{
					localMinQRW = totalMinQRWMO[i];
					minElementIndex = i;
				}
			}
			d = minElementIndex * 0.05;
			cout << "\nМетаоптимизация проведена, оптимальное значение d = " << to_string(d);
			outRWMeta << "\n";
			totalMinQRWMO.clear();
			//////////////////////////////////////////////////////////
			//			И Щ Е М
			cout << "\nПроизводим анализ выборки эвристическими методами:\n";
			for (int i = 0; i < sampleTable.size(); i++)
			{
				localQRandom = n * n;
				localQGreed = n * n;
				localQRW = n * n;
				localQFull = n * n;
				localMinQRandom = n * n;
				localMinQGreed = n * n;
				localMinQRW = n * n;
				localMinQFull = n * n;
				for (int i = 0; i < n2; i++)
				{
					perestanovka.push_back(i);
				}
				vector<int> pRandom(perestanovka), pRW(perestanovka), pFull(perestanovka), pGreed(perestanovka);
				int startQ = compareTables(sampleTable[i], sampleTable2[i], n);
				for (int iteration = 0; iteration < iterCnt; iteration++)
				{
					localQRandom = swapRandomIteration(sampleTable[i], sampleTable2[i], pRandom);
					if (localMinQRandom > localQRandom)
						localMinQRandom = localQRandom;
					localAverageQRandom += localQRandom;
					if (localQRandom == 0)
						localProbabilityRandom += 1.0;
					
					localQRW = swapRandomWeightedIteration(sampleTable[i], sampleTable2[i], pRW, d, startQ);
					if (localMinQRW > localQRW)
						localMinQRW = localQRW;
					localAverageQRW += localQRW;
					if (localQRW == 0)
						localProbabilityRW += 1.0;
					
					localQGreed = swapGreedIteration(sampleTable[i], sampleTable2[i], pGreed, startQ);
					if (localMinQGreed > localQGreed)
						localMinQGreed = localQGreed;
					localAverageQGreed += localQGreed;
					if (localQGreed == 0)
						localProbabilityGreed += 1.0;

					localQFull = swapFullIteration(sampleTable[i], sampleTable2[i], pFull);
					if (localMinQFull > localQFull)
						localMinQFull = localQFull;
					localAverageQFull += localQFull;
					if (localQFull == 0)
						localProbabilityFull += 1.0;
				}
				localAverageQFull /= iterCnt;
				localAverageQRandom /= iterCnt;
				localAverageQRW /= iterCnt;
				localAverageQGreed /= iterCnt;
				localProbabilityFull /= iterCnt;
				localProbabilityGreed /= iterCnt;
				localProbabilityRandom /= iterCnt;
				localProbabilityRW /= iterCnt;

				if (localProbabilityFull > 0.)
				{
					totalProbabilityFull += 1.0;
					QFull.push_back(0);
				}
				else
				{
					QFull.push_back(localMinQFull);
				}
				if (localProbabilityGreed > 0.)
				{
					totalProbabilityGreed += 1.0;
					QGreed.push_back(0);
				}
				else
				{
					QGreed.push_back(localMinQGreed);
				}
				if (localProbabilityRandom > 0.)
				{
					totalProbabilityRandom += 1.0;
					QRandom.push_back(0);
				}
				else
				{
					QRandom.push_back(localMinQRandom);
				}
				if (localProbabilityRW > 0.)
				{
					totalProbabilityRW += 1.0;
					QRW.push_back(0);
				}
				else
				{
					QRW.push_back(localMinQRW);
				}
				localAverageQFull = 0.;
				localAverageQRandom = 0.;
				localAverageQRW = 0.;
				localAverageQGreed = 0.;
				localProbabilityFull = 0.;
				localProbabilityGreed = 0.;
				localProbabilityRandom = 0.;
				localProbabilityRW = 0.;
				perestanovka.clear();
				cout << "X";
			}
			totalProbabilityFull /= sizeOfSample;
			totalProbabilityGreed /= sizeOfSample;
			totalProbabilityRW /= sizeOfSample;
			totalProbabilityRandom /= sizeOfSample;
			for (int i = 0; i < QFull.size(); i++)
			{
				averageQFull += QFull[i];
				averageQRandom += QRandom[i];
				averageQRW += QRW[i];
				averageQGreed += QGreed[i];
			}
			averageQFull /= sizeOfSample;
			averageQGreed /= sizeOfSample;
			averageQRW /= sizeOfSample;
			averageQRandom /= sizeOfSample;
			for (int i = 0; i < QFull.size(); i++)
			{
				dispersionFull += pow((QFull[i] - averageQFull), 2);
				dispersionRandom += pow((QRandom[i] - averageQRandom), 2);
				dispersionRW += pow((QRW[i] - averageQRW), 2);
				dispersionGreed += pow((QGreed[i] - averageQGreed), 2);
			}
			dispersionFull /= sizeOfSample;
			dispersionGreed /= sizeOfSample;
			dispersionRW /= sizeOfSample;
			dispersionRandom /= sizeOfSample;
			sigmaFull = sqrt(dispersionFull);
			sigmaGreed = sqrt(dispersionGreed);
			sigmaRW = sqrt(dispersionRW);
			sigmaRandom = sqrt(dispersionRandom);
			deltaQFull = 1.96 * sigmaFull / sqrt(sizeOfSample);
			deltaQGreed = 1.96 * sigmaGreed / sqrt(sizeOfSample);
			deltaQRW = 1.96 * sigmaRW / sqrt(sizeOfSample);
			deltaQRandom = 1.96 * sigmaRandom / sqrt(sizeOfSample);

			//"total average Q:";
			outRandom  << to_string(averageQRandom) << ";";
			outRW << to_string(averageQRW) << ";";
			outGreed << to_string(averageQGreed) << ";";
			outFull << to_string(averageQFull) << ";";


			//"total probability:";
			outRandom << to_string(totalProbabilityRandom) << ";";
			outRW << to_string(totalProbabilityRW) << ";";
			outGreed << to_string(totalProbabilityGreed) << ";";
			outFull << to_string(totalProbabilityFull) << ";";

			//"Dispersion:";

			outRandom << to_string(dispersionRandom) << ";";
			outRW << to_string(dispersionRW) << ";";
			outGreed << to_string(dispersionGreed) << ";";
			outFull << to_string(dispersionFull) << ";";

			//"Sigma:";
			outRandom << to_string(sigmaRandom) << ";";
			outRW << to_string(sigmaRW) << ";";
			outGreed << to_string(sigmaGreed) << ";";
			outFull << to_string(sigmaFull) << ";";

			//"Delta Q:";
			outRandom << to_string(deltaQRandom) << ";";
			outRW << to_string(deltaQRW) << ";";
			outGreed << to_string(deltaQGreed) << ";";
			outFull << to_string(deltaQFull) << ";";

			//"Interval:";
			outRandom << "[" << to_string(averageQRandom - deltaQRandom) << "; " << to_string(averageQRandom + deltaQRandom) << "]\n";
			outRW << "[" << to_string(averageQRW - deltaQRW) << "; " << to_string(averageQRW + deltaQRW) << "]\n";
			outGreed << "[" << to_string(averageQGreed - deltaQGreed) << "; " << to_string(averageQGreed + deltaQGreed) << "]\n";
			outFull << "[" << to_string(averageQFull - deltaQFull) << "; " << to_string(averageQFull + deltaQFull) << "]\n";

			sampleTable.clear();
			sampleTable2.clear();
			QFull.clear();
			QGreed.clear();
			QRandom.clear();
			QRW.clear();
			totalProbabilityFull = 0;
			totalProbabilityGreed = 0;
			totalProbabilityRW = 0;
			totalProbabilityRandom = 0;
			averageQFull = 0;
			averageQGreed = 0;
			averageQRW = 0;
			averageQRandom = 0;
			dispersionFull = 0;
			dispersionGreed = 0;
			dispersionRW = 0;
			dispersionRandom = 0;
			cout << "\nРасчеты завершены для текущей плотности графов d = " << density;

		}
		for (double d = 0; d < 1.01; d += 0.05)
		{
			outRWMeta << to_string(d) << ",";
		}
		outFull.close();
		outRW.close();
		outRandom.close();
		outGreed.close();
		outRWMeta.close();
	}
	if (type == 2)
	{
		vector<int> answerFull, answerRW, answerRandom, answerGreed;
		ifstream graphs;
		string fileName, tempString2;
		char tempString[100];
		char symbol;
		vector<int> graphString;
		int numberOfGraph = 0;
		vector<vector<int>> graphG, graphG2;
		cout << "\nВыбран второй тип работы - работа с заданными графами.\n";
		cout << "\nВведите название файла(например, graphs.txt) с исходными данными: ";
		cin >> fileName;
		graphs.open(fileName);
		while (!graphs.is_open())
		{
			cout << "Ошибка 2. Введен неверный путь к файлу. Повторите ввод.";
			cin >> fileName;
			graphs.open(fileName);
		}
		while (!graphs.eof())
		{
			graphs.getline(tempString, 100);
			tempString2 = tempString;
			for (int i = 0; i < tempString2.size(); i++)
			{
				if (tempString2.size() > 1)
				{
					if ((tempString2[i] == 'G') && (tempString2[i + 1] == '\''))
					{
						numberOfGraph = 2;
						break;
					}
				}
				else if (tempString2[i] == 'G')
				{
					numberOfGraph = 1;
					break;
				}
				if (numberOfGraph == 0)
				{
					cout << "Ошибка 3. Формат данных в заданном файле отличен от требуемого\n";
					system("pause");
					exit(3);
				}
				if (tempString2[i] == ' ')
					continue;
				else if ((tempString2[i] != '1') && (tempString2[i] != '0'))
				{
					cout << "Ошибка 3. Формат данных в заданном файле отличен от требуемого\n";
					system("pause");
					exit(3);
				}
				else if (tempString2[i] == '1')
					graphString.push_back(1);
				else
					graphString.push_back(0);
			}
			if ((numberOfGraph == 1) && (graphString.size() > 0))
				graphG.push_back(graphString);
			else if ((numberOfGraph == 2) && (graphString.size() > 0))
				graphG2.push_back(graphString);
			graphString.clear();
		}
		printMatrInt(graphG, "G: ");
		printMatrInt(graphG2, "G': ");
		//проверка правильности введенных графов
		//квардатная ли матрица
		if (graphG.size() != graphG[0].size())
		{
			cout << "Ошибка 4. Неверные размерности графов в файле.\n";
			system("pause");
			exit(4);
		}
		else if (graphG.size() != graphG[0].size())
		{
			cout << "Ошибка 4. Неверные размерности графов в файле.\n";
			system("pause");
			exit(4);
		}
		//проверка главной диагонали
		for (int i = 0; i < graphG.size(); i++)
		{
			if (graphG[i][i] != 0)
			{
				cout << "Ошибка 5. Единицы на главной диагонали матрицы графа.\n";
				system("pause");
				exit(5);
			}
		}
		for (int i = 0; i < graphG2.size(); i++)
		{
			if (graphG2[i][i] != 0)
			{
				cout << "Ошибка 5. Единицы на главной диагонали матрицы графа.\n";
				system("pause");
				exit(5);
			}
		}
		//проверка остальных ячеек
		for (int i = 0; i < graphG.size()-1; i++)
		{
			for (int j = i+1; j < graphG[i].size(); j++)
			{
				if (graphG[i][j] != graphG[j][i])
				{
					cout << "Ошибка 6. Значения ячеек матриц смежности не удовлетворяют требованиям.\n";
					system("pause");
					exit(6);
				}
			}
		}
		for (int i = 0; i < graphG2.size()-1; i++)
		{
			for (int j = i+1; j < graphG2[i].size(); j++)
			{
				if (graphG2[i][j] != graphG2[j][i])
				{
					cout << "Ошибка 6. Значения ячеек матриц смежности не удовлетворяют требованиям.\n";
					system("pause");
					exit(6);
				}
			}
		}
		//введенные графы верны
		n2 = graphG.size();
		n = graphG2.size();
		//расчет плотности графа
		double edgesCnt = 0.;
		for (int i = 0; i < graphG.size(); i++)
		{
			for (int j = 0; j < graphG[i].size(); j++)
			{
				if (graphG[i][j] == 1)
				{
					edgesCnt += 1.;
				}
			}
		}
		density = edgesCnt / (n2 * (n2 - 1));
		edgesCnt = 0.;
		cout << "\nДельта = " << (n2 - n) << "\n";
		cout << "\nЧисло различных подграфов при заданных n2 и n: " << (fact(n2) / fact(n) / fact(n2 - n)) << "\n";
		cout << "Плотность графа G: " << density << "\n";
		for (int i = 0; i < graphG2.size(); i++)
		{
			for (int j = 0; j < graphG2[i].size(); j++)
			{
				if (graphG2[i][j] == 1)
				{
					edgesCnt += 1.;
				}
			}
		}
		density = edgesCnt / (n * (n - 1));
		cout << "Плотность графа G': " << density << "\n";
		cout << "Результаты метаоптимизации: "<< d << "\n";
		cout << "Введите число итераций каждого из эвристических алгоритмов(3-999): ";
		cin >> iterCnt;
		while (!((iterCnt > 2) && (iterCnt < 1000)))
		{
			cout << "Ошибка 7. Превышение диапазона данных переменной iterCnt. Повторите ввод.\n";
			cout << "iterCnt = ";
			cin >> iterCnt;
		}
		cout << "\n";
		//////////////////////////////////////////////////////////
		//			И Щ Е М
		cout << "\nПроизводим анализ пары графов эвристическими методами:\n";
		localQRandom = n * n;
		localQGreed = n * n;
		localQRW = n * n;
		localQFull = n * n;
		localMinQRandom = n * n;
		localMinQGreed = n * n;
		localMinQRW = n * n;
		localMinQFull = n * n;
		for (int i = 0; i < n2; i++)
		{
			perestanovka.push_back(i);
		}
		vector<int> pRandom(perestanovka), pRW(perestanovka), pFull(perestanovka), pGreed(perestanovka);
		int startQ = compareTables(graphG, graphG2, n);
		for (int iteration = 0; iteration < iterCnt; iteration++)
		{
			localQRandom = swapRandomIteration(graphG, graphG2, pRandom);
			if (localMinQRandom > localQRandom)
				localMinQRandom = localQRandom;
			localAverageQRandom += localQRandom;
			if (localQRandom == 0)
			{
				localProbabilityRandom += 1.0;
				vector<int> shortPerestanovka;
				for (int i = 0; i < n; i++)
					shortPerestanovka.push_back(pRandom[i]);
				answerRandom = shortPerestanovka;
			}

			localQRW = swapRandomWeightedIteration(graphG, graphG2, pRW, d, startQ);
			if (localMinQRW > localQRW)
				localMinQRW = localQRW;
			localAverageQRW += localQRW;
			if (localQRW == 0)
			{
				localProbabilityRW += 1.0;
				vector<int> shortPerestanovka;
				for (int i = 0; i < n; i++)
					shortPerestanovka.push_back(pRW[i]);
				answerRW = shortPerestanovka;
			}

			localQGreed = swapGreedIteration(graphG, graphG2, pGreed, startQ);
			if (localMinQGreed > localQGreed)
				localMinQGreed = localQGreed;
			localAverageQGreed += localQGreed;
			if (localQGreed == 0)
			{
				localProbabilityGreed += 1.0;
				vector<int> shortPerestanovka;
				for (int i = 0; i < n; i++)
					shortPerestanovka.push_back(pGreed[i]);
				answerGreed = shortPerestanovka;
			}

			localQFull = swapFullIteration(graphG, graphG2, pFull);
			if (localMinQFull > localQFull)
				localMinQFull = localQFull;
			localAverageQFull += localQFull;
			if (localQFull == 0)
			{
				localProbabilityFull += 1.0;
				vector<int> shortPerestanovka;
				for (int i = 0; i < n; i++)
					shortPerestanovka.push_back(pFull[i]);
				answerFull = shortPerestanovka;
			}
		}
		localAverageQFull /= iterCnt;
		localAverageQRandom /= iterCnt;
		localAverageQRW /= iterCnt;
		localAverageQGreed /= iterCnt;
		localProbabilityFull /= iterCnt;
		localProbabilityGreed /= iterCnt;
		localProbabilityRandom /= iterCnt;
		localProbabilityRW /= iterCnt;
		//Среднее качество решения
		cout << "Среднее качество решения за итерацию: \n";
		cout << "\tМетод полного перебора: " << to_string(localAverageQFull) << "\n";
		cout << "\tЖадный метод: " << to_string(localAverageQGreed) << "\n";
		cout << "\tМетод случайного перебора: " << to_string(localAverageQRandom) << "\n";
		cout << "\tМетод случайного взвешенного перебора: " << to_string(localAverageQRW) << "\n";
		
		//Средняя вероятность получения решения
		cout << "\nСредняя вероятность получения решения за итерацию: \n";
		cout << "\tМетод полного перебора: " << to_string(localProbabilityFull) << "\n";
		cout << "\tЖадный метод: " << to_string(localProbabilityGreed) << "\n";
		cout << "\tМетод случайного перебора: " << to_string(localProbabilityRandom) << "\n";
		cout << "\tМетод случайного взвешенного перебора: " << to_string(localProbabilityRW) << "\n";
		cout << "\nРасчеты завершены для текущей плотности графов d = " << density << "\n";

		if (answerFull.size() > 0)
			printVectorInt(answerFull, "Перестановка изоморфного подграфа по методу полного перебора");
		if (answerGreed.size() > 0)
			printVectorInt(answerGreed, "Перестановка изоморфного подграфа по жадному методу");
		if (answerRandom.size() > 0)
			printVectorInt(answerRandom, "Перестановка изоморфного подграфа по методу случайного перебора");
		if (answerRW.size() > 0)
			printVectorInt(answerRW, "Перестановка изоморфного подграфа по методу случайного взвешенного перебора");
	}
	cout << "\n\nВсе расчеты завершены. Результаты расчетов для первого типа работы находятся в файлах в директории с программой";
	cout << "\nВы можете закрыть программу.";
	system("pause");
	return 0;
}