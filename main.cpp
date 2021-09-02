#include <iostream>
#include <utility>
#include <vector>

using namespace std;

#define MAT_SIZE 2
#define ONE 1

// Task1 : Define an affine transformation
class AffineObject
{
private:
    std::vector<std::vector<double>> mat;
    std::vector<std::vector<double>> dia;
public:
    AffineObject(vector<vector<double>> M,vector<vector<double>> d);
    double getm00( ) const;
    double getm01( ) const;
    double getm10( ) const;
    double getm11( ) const;
    double getdx( ) const;
    double getdy( ) const;
};

AffineObject::AffineObject( vector<vector<double>> M, vector<vector<double>> d) : mat(std::move(M)), dia(std::move(d)) {}

double AffineObject::getm00() const
{
    return mat[0][0];
}
double AffineObject::getm01() const
{
    return mat[0][1];
}
double AffineObject::getm10() const
{
    return mat[1][0];
}
double AffineObject::getm11() const
{
    return mat[1][1];
}

double AffineObject::getdx() const
{
    return dia[0][0];
}
double AffineObject::getdy() const
{
    return dia[0][1];
}

class AffineTransformation
{
private:
    double x;
    double y;
    AffineObject tran;
public:
    AffineTransformation(double &&xx, double &&yy, AffineObject t) : x(xx), y(yy), tran(std::move(t)){}
    void transformation();
    void inv_transformation();
    double getx() const;
    double gety() const;
    void print() const;
};

void AffineTransformation::print() const
{
    cout << x << " " << y << endl;
}

double AffineTransformation::getx() const
{
    return x;
}

double AffineTransformation::gety() const
{
    return y;
}

// Task2 : Implement a transformation function
void AffineTransformation::transformation()
{
    double tx = x;
    double ty = y;
    x = tx * tran.getm00() + ty * tran.getm01() + tran.getdx();
    y = tx * tran.getm10() + ty * tran.getm11() + tran.getdy();
}

// Task3 : Inverse transformations
void AffineTransformation::inv_transformation()
{
    double m00,m01,m10,m11,dx,dy;
    m00 = tran.getm00();
    m01 = tran.getm01();
    m10 = tran.getm10();
    m11 = tran.getm11();
    dx = tran.getdx();
    dy = tran.getdy();

    double denom = m00 * m11 - m01 * m10;

    if (denom != 0)
    {
        // calculate the inverse values
        double im00,im01,im10,im11,idx,idy;
        im00 = m11 / denom;
        im01 = -m01 / denom;
        im10 = -m10 / denom;
        im11 = m00 / denom;
        idx = - dx * im00 - dy * im01;
        idy = - dx * im10 - dy * im11;

        double tx = x;
        double ty = y;

        x = tx * im00 + ty * im01 + idx;
        y = tx * im10 + ty * im11 + idy;
    }
    else
    {
        cout << "Cannot calculate the inverse" << endl;
    }
};

class AffineCTransformation
{
private:
    AffineObject a;
    AffineObject b;
public:
    AffineCTransformation(AffineObject aa, AffineObject bb) : a(std::move(aa)), b(std::move(bb)) {};
    AffineObject comb_transformation();
};

// Task4 : Combine transformations
AffineObject AffineCTransformation::comb_transformation()
{

    std::vector<std::vector<double>> cM(MAT_SIZE,std::vector<double>(MAT_SIZE,0));
    std::vector<std::vector<double>> cd(ONE,std::vector<double>(MAT_SIZE,0));

    cM[0][0] = a.getm00() * b.getm00() + a.getm01() * b.getm10();
    cM[0][1] = a.getm00() * b.getm01() + a.getm01() * b.getm11();
    cM[1][0] = a.getm10() * b.getm00() + a.getm11() * b.getm10();
    cM[1][1] = a.getm10() * b.getm01() + a.getm11() * b.getm11();
    cd[0][0] = a.getdx() + a.getm00() * b.getdx() + a.getm01() * b.getdy();
    cd[0][1] = a.getdy() + a.getm10() * b.getdx() + a.getm11() * b.getdy();

    AffineObject ctran(cM, cd);
    return ctran;
};

// example
int main()
{
    std::vector<std::vector<double>> M1(MAT_SIZE,std::vector<double>(MAT_SIZE,0));
    std::vector<std::vector<double>> M2(MAT_SIZE,std::vector<double>(MAT_SIZE,0));

    for(int i = 0; i < MAT_SIZE; i++)
    {
        for(int j = 0; j < MAT_SIZE; j++)
        {
            std::cin >> M1[i][j];
        }
    }

    for(int i = 0; i < MAT_SIZE; i++)
    {
        for(int j = 0; j < MAT_SIZE; j++)
        {
            std::cin >> M2[i][j];
        }
    }

    std::vector<std::vector<double>> d1(MAT_SIZE,std::vector<double>(MAT_SIZE,0));
    std::vector<std::vector<double>> d2(MAT_SIZE,std::vector<double>(MAT_SIZE,0));

    for(int i = 0; i < ONE; i++)
    {
        for(int j = 0; j < MAT_SIZE; j++)
        {
            std::cin >> d1[i][j];
        }
    }

    for(int i = 0; i < ONE; i++)
    {
        for(int j = 0; j < MAT_SIZE; j++)
        {
            std::cin >> d2[i][j];
        }
    }

    AffineObject examp(M1, d1);
    AffineObject examp1(M2, d2);

    AffineTransformation AF(1,0,examp);

    cout << "Initial values:" << endl;
    AF.print();

    AF.transformation();
    cout << "Transformation:" << endl;
    AF.print();

    AF.inv_transformation();
    cout << "Inverse Transformation:" << endl;
    AF.print();

    AffineCTransformation ACF(examp,examp1);
    AffineObject comb = ACF.comb_transformation();

    AffineTransformation AF2(AF.getx(),AF.gety(),comb);
    cout << "Combined Transformation:" << endl;
    AF.print();

    return 0;
}
