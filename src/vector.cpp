#include "vector.hpp"
#include "geometry.hpp"

using namespace std;

/* Returns a vector vzip that is the iteration of two other vectors v1 and v2, both with size n=n1=n2,
    / where the first element of vzip is the iteration of the first elements of v1 and v2 together,
    / the second element of vzip is the iteration of the second elements of v1 anv v2 together, and so on.
    / returns: a vector with n vectors */
std::vector<std::vector<double> > zip_vector(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<std::vector<double> > vzip;

    if (v1.size() == v2.size())
    {
        for (int i = 0; i < v1.size(); i++)
        {
            vzip.insert(vzip.end(), {v1[i], v2[i]});
        }
    }
    else
    {
        cout << "Error: Vectors have different sizes and cannot be zipped together." << '\n';
    }

    return vzip;
}

/* Same as previous function, but for a vector of vectors of vectors */
std::vector<std::vector<std::vector<double> > > zip_matrix(std::vector<std::vector<double> > vertices1, std::vector<std::vector<double> > vertices2)
{
    int n1 = vertices1.size();
    int n2 = vertices2.size();
    std::vector<std::vector<std::vector<double> > > mzip;

    if (n1 == n2)
    {
        for (int i = 0; i < n1; i++)
        {
            mzip.insert(mzip.end(), {vertices1[i], vertices2[i]});
        }
    }
    else
    {
        cout << "Error: Vertex matrices have different sizes and cannot be zipped together." << '\n';
    }

    return mzip;
}

// Scales vector by a constant
std::vector<double> scale_vector(std::vector<double> vector, double s)
{
    std::vector<double> scaled_vector;
    for (int i = 0; i < vector.size(); i++)
    {
        scaled_vector.push_back(vector[i] * s);
    }
    return scaled_vector;
}

// Summation of vectors: creates a vector containing the sums of the elements of two other vectors
std::vector<double> add_vectors(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size() == v2.size());

    std::vector<double> vector_summation;
    vector_summation.reserve(v1.size());

    std::transform(v1.begin(), v1.end(), v2.begin(),
                   std::back_inserter(vector_summation), std::plus<double>());

    return vector_summation;
}

// Summation of vectors: creates a vector containing the sums of the elements of four other vectors
std::vector<double> add_vectors(const std::vector<double> &v1, const std::vector<double> &v2, const std::vector<double> &v3, const std::vector<double> &v4)
{
    assert(v1.size() == v2.size());
    assert(v3.size() == v4.size());
    assert(v4.size() == v1.size());

    std::vector<double> vector_summation;
    vector_summation.reserve(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        vector_summation.push_back(v1[i] + v2[i] + v3[i] + v4[i]);
    }

    return vector_summation;
}

// Summation of vectors: creates a vector containing the sums of the elements of two other vectors
std::vector<double> subtract_vectors(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size() == v2.size());

    std::vector<double> vector_subtraction;
    vector_subtraction.reserve(v1.size());

    std::transform(v1.begin(), v1.end(), v2.begin(),
                   std::back_inserter(vector_subtraction), std::minus<double>());

    return vector_subtraction;
}

// Calculates the dot product between two vectors v1 and v2
double get_dot_product(std::vector<double> v1, std::vector<double> v2)
{
    double dot_product = 0;

    for (int i = 0; i < v1.size(); i++)
    {
        dot_product += v1[i] * v2[i];
    }
    return dot_product;
}

// Calculates the dot product between a matrix m1 (2x2) and a vector v1 (2D)
std::vector<double> get_dot_product_matrix(std::vector<std::vector<double> > m1, std::vector<double> v1)
{
    std::vector<double> dot_product;

    for (int i = 0; i < m1.size(); i++)
    {
        double product = 0;
        for (int j = 0; j < v1.size(); j++)
        {
            product += m1[i][j] * v1[j];
        }
        dot_product.push_back(product);
    }
    return dot_product;
}

// Calculates the unit vector between two vectors
std::vector<double> get_unit_vector(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> vector;
    vector.insert(vector.end(), {(v1[0] - v2[0]), (v1[1] - v2[1])});

    double dist = get_euclidian_distance(v1[0], v1[1], v2[0], v2[1]);
    std::vector<double> uv;
    uv.insert(uv.end(), {(vector[0] / dist), (vector[1] / dist)});

    return uv;
}

// Converts angle theta to unit vector
std::vector<double> angle_2_vector(double theta)
{
    double x = cos(theta);
    double y = sin(theta);

    // Convert to unit vector
    std::vector<double> v1{x, y};
    std::vector<double> v2{0, 0};

    std::vector<double> uv = get_unit_vector(v1, v2);

    return uv;
}

// Converts vector to angle
double vector_2_angle(double x, double y)
{
    return atan2(y, x);
}

// Calculates the sqrt of the square of each element in the vector
std::vector<double> absolute_vector(std::vector<double> vector)
{
    std::vector<double> absolute_vector;
    for (int i = 0; i < vector.size(); i++)
    {
        absolute_vector.push_back(fabs(vector[i]));
    }
    return absolute_vector;
}

/*// Calculates the cross product between two vectors v1 and v2
std::vector<double> get_cross_product(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> cross_product{(v1[1] * v2[2] - v1[2] * v2[1]),
                                      (v1[2] * v2[0] - v1[0] * v2[2]),
                                      (v1[0] * v2[1] - v1[1] * v2[0])};

    return cross_product;
}*/