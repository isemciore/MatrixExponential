#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cmath>

typedef double matrix_VT;
typedef std::vector<matrix_VT> rowType;
class Matrix
{
 public:
    typedef unsigned int index;
    Matrix( );
    explicit Matrix(const std::size_t &, const std::size_t &);
    Matrix( const Matrix& ); //copy
    explicit Matrix(const std::size_t &);
    ~Matrix( );
    
    Matrix& operator= ( const Matrix& ); //Assigncopy
    Matrix operator+ ( const Matrix& ) const;
    Matrix operator* ( const Matrix& ) const; //Mtrx mult
    Matrix operator* ( const matrix_VT & ) const;
    Matrix operator-( const Matrix& ) const;
    Matrix operator-( ) const;   
    friend Matrix operator* (const matrix_VT & , const Matrix& );
    
    Matrix& operator+=(const Matrix&);
    Matrix& operator*=(const Matrix&);
    double norm() const;
    Matrix& transpose( );
    
    rowType& operator[](const index & i );
    const rowType operator[](const index & i ) const;
	double maxElement() const;
    std::size_t rows() const;
    std::size_t cols() const;
    void set_id();
    friend std::ostream& operator<<(std::ostream& os, const Matrix&);
	friend Matrix expMatrix(Matrix&,std::size_t);
	void printMatrix(std::string basic_string);

protected:
 private:
    std::vector<std::vector<matrix_VT>>       m_vectors;
    std::size_t                 m_rows;
    std::size_t                 m_cols;
    double						tol = 5e-10;

    friend std::istream& operator>> ( std::istream&, Matrix& );
};

Matrix expMatrix(Matrix& src, std::size_t tot_it = 1000){
    if( src.cols() != src.rows()){
        throw std::out_of_range("not square matrix");
    }
    double maxElt = src.maxElement();
    double tempMax = maxElt;

    Matrix temp(src.rows()); //saves the current sequence
    temp.set_id();

    Matrix A = src;
    Matrix answer = temp;//Output, temp is the first term in sum
    double constant = 1;

    for(std::size_t i = 1; i < tot_it;i++){
        constant = i;
        temp *= (1/constant)*A; //Calculates the next term (A^n/n!)
        answer += temp;//f(x) = f(x) + next term

        //Checking tolerance/stopping criteria
        double decr = maxElt/i; //once less then we know the sequence begins decrease
        tempMax *= decr;
        if( (decr < 0.1) && (tempMax < src.tol)){
            break;
        }
        if(i == 999){
            std::cout << "Reached 1000 iteration, stopping \n";
        }
    }

    return answer;
}
std::istream& operator>> ( std::istream& inStream, Matrix& src){
	src.m_vectors.clear();
	std::string temp;
	inStream >> temp;
	std::size_t rowNo = 0;
	std::size_t colNo = 0;
	if(temp  =="["){
		rowType dummy;
		src.m_vectors.push_back(dummy);
		while(inStream>>temp){
			if(temp == "]"){
				break;
			}
			else if(temp == ";"){
				colNo = 0;
				rowNo++;
				rowType dummy;
				src.m_vectors.push_back(dummy);
			}
			else{
				std::cout << temp << "\n";
				src.m_vectors[rowNo].push_back(std::stoi(temp));
				colNo++;
			}
		}
	}
	src.m_cols = colNo;
	src.m_rows = rowNo +1;
	return inStream;
}

std::ostream& operator<< (std::ostream& output, const Matrix & src){
	output << "[ ";
	for(std::size_t rowNum = 0; rowNum < src.m_rows; rowNum++){
		for(std::size_t colNum = 0; colNum < src.m_cols; colNum++){
			output << src.m_vectors[rowNum][colNum]<<" ";
		}
		if(rowNum < src.m_rows-1){
			output << "\n; ";
		}
		else{
			output << "]";	
		}
	}
    output << "\n";
	return output;
}




void Matrix::set_id(){
    if(m_cols!=m_rows){
        throw std::out_of_range("not square");
    }
    for(std::size_t i = 0; i < m_rows;i++){
        for(std::size_t j = 0; j<m_cols;j++){
            if(i==j) {
                m_vectors[i][j] = 1;
            }
            else{
                m_vectors[i][j] = 0;
            }
        }
    }
}
//Matrix operator* ( int, const Matrix& ); //??

//construct rektangular mtrx
Matrix::Matrix(const std::size_t & height, const std::size_t & width)
	: m_vectors(height,std::vector<matrix_VT>(width))
	, m_rows(height)
	, m_cols(width)
	{
}


Matrix::Matrix(const std::size_t & sideLength)
	: m_vectors(sideLength,std::vector<matrix_VT>(sideLength))
	, m_rows(sideLength)
	, m_cols(sideLength)
	{
}
//Copy construct
Matrix::Matrix(const Matrix& src)
	: m_vectors(src.m_vectors)
	, m_rows(src.m_rows)
	, m_cols(src.m_cols)
	{
}


//Copy via assignment 
Matrix& Matrix::operator= ( const Matrix& src){
	if(this==&src) return *this;
	
	m_rows = src.m_rows;
	m_cols = src.m_cols;
	m_vectors = src.m_vectors; //Assuming vector has copy assignable
	return *this;	
}




Matrix& Matrix::transpose( ){
	
	Matrix tmp_mtrx(m_cols,m_rows);
	//int dummy = 0;
	for(std::size_t i = 0; i < m_rows; i++){
		for(std::size_t j = 0; j < m_cols; j++){
			//dummy += m_vectors[i][j];
			//std::cout << m_rows << " " <<m_cols<<" *.'\n";
			//std::cout << i << " " << j <<" *.*\n";
			tmp_mtrx.m_vectors[j][i] = m_vectors[i][j];
		}
	}
	m_vectors = tmp_mtrx.m_vectors;
	//std::cout << dummy << "  \n";
	std::size_t tmpSize = m_rows;
	m_rows = m_cols;
	m_cols = tmpSize;
	
	return *this;	
}




void Matrix::printMatrix(std::string fileName) {
	std::ofstream myfile;
    myfile.precision(15);
	myfile.open(fileName,std::fstream::out);
	if(myfile.is_open()){
		for(std::size_t rowNum =0; rowNum < m_rows; rowNum++){
			for(std::size_t colNum = 0; colNum < m_cols; colNum++){
				myfile << m_vectors[rowNum][colNum] << " ";
			}
			myfile << "\n";
		}
		myfile << "\n";
		myfile.close();
	} else{
		std::cout << "Failed to open file";
	}
}








Matrix Matrix::operator+ (const Matrix& addMe) const{ //Addition 
	if((m_rows != addMe.m_rows)|| (m_cols != addMe.m_rows)){
		throw std::out_of_range("matrices size does not match");
	}
	
	Matrix this_obj = *this;//
	for(std::size_t i = 0; i < m_rows; i++){
		for(std::size_t j = 0; j < m_cols; j++){
			this_obj.m_vectors[i][j] += addMe.m_vectors[i][j];	
		}	
	}
	return this_obj;
}

//kopierar sig sj. ändrar bara på kopian, aldrig sig själv
Matrix Matrix::operator* ( const matrix_VT & k) const { //sclar mult
	Matrix this_obj = *this;
	for(std::size_t i = 0; i < m_rows; i++){
		for(std::size_t j = 0; j < m_cols; j++){
			this_obj.m_vectors[i][j] *= k;	
		}	
	}
	return this_obj;
}

Matrix Matrix::operator* ( const Matrix& multMe) const{ //mtrx mult
	if((m_rows != multMe.m_cols)|| (m_cols != multMe.m_rows)){
		throw std::out_of_range("matrices size does not match");
	}
	
	Matrix temp(m_rows,multMe.m_cols);
	matrix_VT tmp_sum = 0;
	for(std::size_t i=0; i<m_rows; i++){
		for(std::size_t j=0; j < m_cols; j++){
			tmp_sum = 0;
			for(std::size_t itNum = 0; itNum < m_rows; itNum++){
				tmp_sum += m_vectors[i][itNum]*multMe.m_vectors[itNum][j];
			}
			temp[i][j] = tmp_sum;
		}
	}	
	return temp;
}

Matrix operator* (const matrix_VT & scalar ,const Matrix& matrix){
	std::size_t nRows = matrix.rows();
	std::size_t nCols = matrix.cols();
	Matrix temp(nRows, nCols);
	for(std::size_t i = 0; i < nRows;i++){
		for(std::size_t j = 0; j < nCols; j++){	
			temp[i][j] = scalar*matrix[i][j];
		}
	}
	
	return temp;	
}
    
    
    
    
Matrix& Matrix::operator+=(const Matrix& addMe){
	if((m_cols != addMe.m_cols)|| (m_rows != addMe.m_rows)){
		throw std::out_of_range("matrices size does not match");
	}
	for(std::size_t i = 0; i < m_rows;i++){
		for(std::size_t j = 0; j < m_cols; j++){	
			m_vectors[i][j] += addMe[i][j];
		}
	}
	return *this;
}

Matrix& Matrix::operator*=(const Matrix& multThis){
	if((m_rows != multThis.m_cols)|| (m_cols != multThis.m_rows)){
		throw std::out_of_range("matrices size does not match");
	}
	Matrix temp(m_rows,m_cols);
	matrix_VT tmp_sum = 0;
	for(std::size_t i=0; i<m_rows; i++){
		for(std::size_t j=0; j < m_cols; j++){
			tmp_sum = 0;
			for(std::size_t itNum = 0; itNum < m_rows; itNum++){
				tmp_sum += m_vectors[i][itNum]*multThis[itNum][j];
			}
			temp[i][j] = tmp_sum;
		}
	}
	for(std::size_t i=0; i < m_rows;i++){
		for(std::size_t j=0;j<m_cols;j++){
			m_vectors[i][j] = temp[i][j];
		}
	}
	return * this;
}

double Matrix::norm()const{
	//L2,1 norm
	double tmpAllsum = 0;
	for(std::size_t i=0; i <m_rows; i++){
		double tmpRowSum = 0;
		for(std::size_t j=0; j<m_cols;j++){
			tmpRowSum += m_vectors[i][j]*m_vectors[i][j];
		}
		tmpAllsum += sqrt(tmpRowSum);
	}
	return tmpAllsum;
}


double Matrix::maxElement() const {
    double maxValue = 1;
    for(std::size_t i = 0; i <m_rows;i++){
        for(std::size_t j = 0; j < m_cols; j++){
            if( fabs(m_vectors[i][j]) > maxValue){
                maxValue = m_vectors[i][j];
            }
        }
    }
    return  fabs(maxValue);
}
    

Matrix Matrix::operator-( const Matrix& subtrMe) const{
	if((m_rows != subtrMe.m_rows)|| (m_cols != subtrMe.m_rows)){
		throw std::out_of_range("matrices size does not match");
	}
	Matrix this_obj = *this;
	for(std::size_t i = 0; i < m_rows; i++){
		for(std::size_t j = 0; j < m_cols; j++){
			this_obj.m_vectors[i][j] -= subtrMe.m_vectors[i][j];	
		}	
	}
	return this_obj;
	
}


Matrix Matrix::operator-( ) const{
	Matrix copy_of_this = *this;
	for(std::size_t i = 0; i < m_rows; i++){
		for(std::size_t j = 0; j < m_cols; j++){
			copy_of_this.m_vectors[i][j] = -m_vectors[i][j];
		}
	}
	return copy_of_this;	
}

Matrix::~Matrix(){
}

rowType& Matrix::operator[](const index & i ){

	return m_vectors[i];
}

const rowType Matrix::operator[](const index & i ) const{
	const rowType tmp = m_vectors[i];
	return tmp;
}

std::size_t Matrix::rows() const{
	return m_rows;
}
std::size_t Matrix::cols() const{
	return m_cols;
}


#endif // MATRIX_H

