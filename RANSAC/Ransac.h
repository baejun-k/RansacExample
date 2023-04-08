#ifndef __RANSAC_HEADER__
#define __RANSAC_HEADER__


#include <vector>
#include <array>
#include <random>


template<typename ElementType>
struct Point2;

template<>
struct Point2<float>
{
	float x;
	float y;
};

template<>
struct Point2<double>
{
	double x;
	double y;
};

template<typename ElementType>
class Ransac2D
{
public:
	using Point = Point2<ElementType>;

private:
	std::array<ElementType, 2> EstimateModel(const std::vector<const Point*>& dataPtrs) const
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<size_t> distrib(0, dataPtrs.size() - 1);

		const size_t maxNumData = dataPtrs.size();
		const Point2<ElementType>& p1 = *dataPtrs[distrib(gen)];
		const Point2<ElementType>& p2 = *dataPtrs[distrib(gen)];

		ElementType a = (p2.y - p1.y) / (p2.x - p1.x);
		ElementType b = p1.y - a * p1.x;

		return { a, b };
	}

	ElementType MeasurePerformance(const std::vector<const Point*>& dataPtrs,
								   const std::array<ElementType, 2>& model)
	{
		ElementType error = { 0 };
		const ElementType sqrt = std::sqrt(model[0] * model[0] + 1);
		for (const Point* p : dataPtrs)
		{
			ElementType dist = std::abs(p->y - model[0] * p->x - model[1]) / sqrt;
			error += dist;
		}
		return (error / dataPtrs.size());
	}

	void Suffle(std::vector<const Point*>& dataPtrs) const
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::shuffle(dataPtrs.begin(), dataPtrs.end(), gen);
	}

public:
	std::array<ElementType, 2> Calc(const std::vector<Point2<ElementType>>& data,
									const size_t& maxIterations,
									const ElementType& threshold,
									const size_t& sampleSize)
	{
		std::vector<const Point*> dataPtrs(data.size());
		for (size_t i = 0; i < data.size(); ++i)
		{
			dataPtrs[i] = &data[i];
		}

		std::array<ElementType, 2> bestModel = { 0 };
		ElementType bestPerformance = std::numeric_limits<ElementType>::max();

		std::vector<const Point*> sample;
		std::vector<const Point*> inliers;
		std::vector<const Point*> outliers;
		sample.reserve(data.size());
		inliers.reserve(data.size());
		outliers.reserve(data.size());

		for (size_t i = 0; i < maxIterations; ++i)
		{
			sample.clear();
			inliers.clear();
			outliers.clear();

			Suffle(dataPtrs);

			for (size_t k = 0; k < sampleSize; ++k)
			{
				sample.push_back(dataPtrs[k]);
			}

			const std::array<ElementType, 2> model = EstimateModel(sample);
			const ElementType modelSqrt = std::sqrt(model[0] * model[0] + 1);

			// 이상치 제외한 데이터를 사용하여 모델 재추정
			for (const Point* p : dataPtrs)
			{
				ElementType dist = std::abs(p->y - model[0] * p->x - model[1]) / modelSqrt;
				if (threshold < dist)
				{
					outliers.push_back(p);
				}
				else
				{
					inliers.push_back(p);
				}
			}
			const std::array<ElementType, 2> newModel = EstimateModel(inliers);

			ElementType performance = MeasurePerformance(dataPtrs, newModel);
			if (performance < bestPerformance)
			{
				bestModel = newModel;
				bestPerformance = performance;
			}
		}

		return bestModel;
	}

};


#endif // !__RANSAC_HEADER__