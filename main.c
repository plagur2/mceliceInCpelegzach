#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// Function to open the "Choose File" dialog
#define ROWS 4
#define COLS 8
#define ROWS_G 4
#define COLS_G 8
#define ROWS_H 4
#define COLS_H 8

typedef struct PrivateKey {
    int(*G)[8];
    int(*H)[8];
    int(*S)[4];
    int(*P)[8];
    int(*GPRIME)[8];
}privateKey1,*privateKEYPTR;



int length=0;
int n_G;

int p_G;

int col_G;
#include <stdio.h>

#define N 8

// Function to get the cofactor of an element
void getCofactor(int(*mat)[8], int temp[N][N], int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// Recursive function to find the determinant
int determinantOfMatrix(int(*mat)[8], int n) {
    int D = 0;
    if (n == 1)
        return mat[0][0];

    int temp[N][N];
    int sign = 1;

    for (int f = 0; f < n; f++) {
        getCofactor(mat, temp, 0, f, n);
        D += sign * mat[0][f] * determinantOfMatrix(temp, n - 1);
        sign = -sign;
    }

    return D;
}

void getCofactor4(int(*mat)[4], int temp[4][4], int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// Recursive function to find the determinant
int determinantOfMatrix4(int(*mat)[4], int n) {
    int D = 0;
    if (n == 1)
        return mat[0][0];

    int temp[4][4];
    int sign = 1;

    for (int f = 0; f < n; f++) {
        getCofactor4(mat, temp, 0, f, n);
        D += sign * mat[0][f] * determinantOfMatrix4(temp, n - 1);
        sign = -sign;
    }

    return D;
}

// Function to find the adjoint of the matrix
void adjoint(int(*mat)[8], int adj[8][8]) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            int temp[8][8];
            getCofactor(mat, temp, i, j, 8);
            adj[j][i] = determinantOfMatrix(temp, 8 - 1); // Transpose indices
            adj[j][i] *= (i + j) % 2 == 0 ? 1 : -1; // Apply sign
        }
    }
}

// Function to find the inverse of the matrix (rounded to integers)
void inverse8(int(*mat)[8], int(*inv)[8]) {
    int det = determinantOfMatrix(mat, 8);
    if (det == 0) {
        printf("Matrix is singular; inverse does not exist.\n");
        return;
    }

    int adj[8][8];
    adjoint(mat, adj);

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            inv[i][j] = adj[i][j] / det; // Integer division
        }
    }
}


void freeMatrix(int** m, int row)
{
    for(int i=0;i<row;i++)
    {
        free(m[i]);
    }
    free(m);
}
int** init( int row,int col)
{
    int (*A)[col] = malloc(sizeof(int[row][col]));
    for (int i = 0; i <4 ; ++i) {
        for (int j = 0; j < 8; ++j) {
            A[i][j] = 0;
            printf("A in %d,%d is: %d \n",i,j,A[i][j]);
        }
    }
    // int** a  = A;
    //  return a;
}
void multiplyMatrix(int(*result)[p_G], int(*m1)[n_G], int(*m2)[p_G], int m, int n, int p) {
    // Initialize all elements to zero
    // Allocate memory for the 2D array

    if(!result){
        return;
    }

    // Initialize all elements to zero
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            result[i][j] = 0;
        }
    }


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {

                result[i][j] += m1[i][k] * m2[k][j];

            }

        }

    }


}
void bitFlip(int *C, int size, int n) {
    if (n == 0) {
        return; // No flip
    }

    if (n == -1) {
        srand(time(NULL)); // Initialize random seed
        n = rand() % size; // Random index between 0 and size - 1
    }

    if (C[n - 1] == 1) {
        C[n - 1] = 0;
    } else {
        C[n - 1] = 1;
    }
}

void genSMatrix(privateKEYPTR p1,int k) {


    while (1) {
        // Generate a random matrix with binary entries (0 or 1)
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                p1->S[i][j] = rand()%2;
            }
        }

        // Check if the matrix is invertible
        int det = determinantOfMatrix4(p1->S, 4);
        if (det != 0.0) {
            return; // Matrix is invertible
        }
    }
}
void genPMatrix(privateKEYPTR p1,int n, int keep) {


    // Initialize the identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            p1->P[i][j] = (i == j) ? 1 : 0;
        }
    }

    // Apply a random permutation (Fisher-Yates shuffle)
    if (!keep) {
        for (int i = n - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            if (i != j) {
                for (int k = 0; k < n; k++) {
                    int temp = p1->P[i][k];
                    p1->P[i][k] = p1->P[j][k];
                    p1->P[j][k] = temp;
                }
            }
        }
    }

}
void modTwoMatrix(int(*matrix)[col_G],int row,int col) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            matrix[i][j] %= 2;
        }
    }
}
int chooseFile(char* filePath) {
    OPENFILENAME ofn; // Structure for file dialog
    char szFile[MAX_PATH]; // Buffer to store the selected file path

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL; // Handle to the parent window (or NULL)
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "All Files\0*.*\0"; // Filter for file types
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn)) {
        strcpy(filePath, szFile); // Store the selected file path
        return 1; // File chosen successfully
    }
    return 0; // User canceled the dialog
}

// Function to read the file and convert it to an array of bits
int* readFileAndConvertToBits(const char* filePath) {
    FILE* file = fopen(filePath, "rb");
    if (!file) {
        printf("Error opening file.\n");

    }
    int* arr =(int*) malloc(4);
    // Read the file byte by byte
    int index=0;
    unsigned char byte;
    while (fread(&byte, sizeof(byte), 1, file) == 1) {
        // Process each bit in the byte
        for (int i = 0; i < 8; i++) {
            int bit = (byte >> i) & 1;
            arr[index] = bit;
            int* arrTEMP = (int*)realloc(arr, (index + 2) * sizeof(int));
            if(arrTEMP)
            {
                arr = arrTEMP;
            }
            index++;
            // Store the bit in your array or perform any other operation
        }
    }
    length = index;
    fclose(file);
    return arr;

}
void gprime(struct PrivateKey* privateKey)
{
    int (*gprime)[8] = malloc(sizeof(int[4][8]));;
    p_G = 8;
    n_G = 4;
    multiplyMatrix(gprime,privateKey->S,  privateKey->G,4,4,8);
    multiplyMatrix(privateKey->GPRIME,gprime,privateKey->P,4,8,8);
    free(gprime);
}
int all_zeros(int* d, int size) {
    for (int i = 0; i < size; ++i) {
        if (d[i] != 0) {
            return 0; // Not all zeros
        }
    }
    return 1; // All zeros
}

// Function to perform syndrome lookup
int syndromeLookup(int(*H)[8], int(*d)[8]) {
    int t[8][4]; // Assuming H is a 4x8 matrix
    int s[24];

    // Convert d to transpose
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 4; ++j) {
            t[i][j] = d[j][i];
        }
    }

    int sCount = 0;
    printf("\n s syndrome \n");
    // Convert d to a 1D array
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 4; ++j) {
            s[sCount] = t[i][j];
            printf("%d,", s[sCount]);
            sCount++;
        }

    }
    printf("\n s syndrome \n");

    if (all_zeros(s, 24)) {
        return 0;
    }

    for (int i = 0; i < 8; ++i) {
        int match = 1;
        for (int j = 0; j < 4; ++j) {
            if (s[j] != t[i][j]) {
                match = 0;
                break;
            }
        }
        if (match) {
            return i + 1;
        }
    }

    return 0; // Not found in syndrome table
}
void InvertMatrix(const double m[16], int invOut[16]) {
    double inv[16], det;
    int i;

    inv[0] = m[5] * m[10] * m[15] -
             m[5] * m[11] * m[14] -
             m[9] * m[6] * m[15] +
             m[9] * m[7] * m[14] +
             m[13] * m[6] * m[11] -
             m[13] * m[7] * m[10];

    inv[4] = -m[4] * m[10] * m[15] +
             m[4] * m[11] * m[14] +
             m[8] * m[6] * m[15] -
             m[8] * m[7] * m[14] -
             m[12] * m[6] * m[11] +
             m[12] * m[7] * m[10];

    inv[8] = m[4] * m[9] * m[15] -
             m[4] * m[11] * m[13] -
             m[8] * m[5] * m[15] +
             m[8] * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4] * m[9] * m[14] +
              m[4] * m[10] * m[13] +
              m[8] * m[5] * m[14] -
              m[8] * m[6] * m[13] -
              m[12] * m[5] * m[10] +
              m[12] * m[6] * m[9];

    inv[1] = -m[1] * m[10] * m[15] +
             m[1] * m[11] * m[14] +
             m[9] * m[2] * m[15] -
             m[9] * m[3] * m[14] -
             m[13] * m[2] * m[11] +
             m[13] * m[3] * m[10];

    inv[5] = m[0] * m[10] * m[15] -
             m[0] * m[11] * m[14] -
             m[8] * m[2] * m[15] +
             m[8] * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0] * m[9] * m[15] +
             m[0] * m[11] * m[13] +
             m[8] * m[1] * m[15] -
             m[8] * m[3] * m[13] -
             m[12] * m[1] * m[11] +
             m[12] * m[3] * m[9];

    inv[13] = m[0] * m[9] * m[14] -
              m[0] * m[10] * m[13] -
              m[8] * m[1] * m[14] +
              m[8] * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1] * m[6] * m[15] -
             m[1] * m[7] * m[14] -
             m[5] * m[2] * m[15] +
             m[5] * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0] * m[6] * m[15] +
             m[0] * m[7] * m[14] +
             m[4] * m[2] * m[15] -
             m[4] * m[3] * m[14] -
             m[12] * m[2] * m[7] +
             m[12] * m[3] * m[6];

    inv[10] = m[0] * m[5] * m[15] -
              m[0] * m[7] * m[13] -
              m[4] * m[1] * m[15] +
              m[4] * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0] * m[5] * m[14] +
              m[0] * m[6] * m[13] +
              m[4] * m[1] * m[14] -
              m[4] * m[2] * m[13] -
              m[12] * m[1] * m[6] +
              m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
             m[1] * m[7] * m[10] +
             m[5] * m[2] * m[11] -
             m[5] * m[3] * m[10] -
             m[9] * m[2] * m[7] +
             m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
              m[0] * m[7] * m[9] +
              m[4] * m[1] * m[11] -
              m[4] * m[3] * m[9] -
              m[8] * m[1] * m[7] +
              m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return;

    det = 1.0 / det;

    for (i = 0; i < 16; i++) {
        double temp = inv[i] * det;
        if (temp == -0.000000){
            invOut[i] = 0;
        }
        else if(temp == -1.000000)
        {
            invOut[i] = -1;
        }
        else if(temp == 1.000000)
        {
            invOut[i] = 1;
        }
    }
    return;
}

void inverse(int(*mat)[4],int(*inv)[4])
{
    double vec[16];
    int inv1[16];
    int count = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j <4 ; ++j) {
            vec[count] = mat[i][j];
            count++;
        }
    }
    count = 0;
    InvertMatrix(vec, inv1);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            inv[i][j] = inv1[count];
            count++;
        }
    }
}

void decrypt( privateKEYPTR p1,int(*c)[8],int(*res)[4]) {
    int (*cHat)[8] = malloc(sizeof(int[1][8]));
    int(*inversedP)[8] = malloc(sizeof(int[8][8]));
    inverse8(p1->P, inversedP);
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            inversedP[i][j] = abs(inversedP[i][j]);
        }
    }
    modTwoMatrix(inversedP, 8, 8);
    int (*transpose)[1] = malloc(sizeof(int[8][1]));

// Compute the transpose
    n_G = 8;
    p_G = 1;
    multiplyMatrix(cHat, c, inversedP, 1, 8, 8);

    for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < 8; ++j) {
            transpose[j][i] = cHat[i][j];
        }
    }
    int (*CH)[1] = malloc(sizeof(int[4][1]));
    n_G = 8;
    p_G = 1;
    multiplyMatrix(CH, p1->H, transpose, 4, 8, 1);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 1; ++j) {
            CH[i][j] = abs(CH[i][j] % 2);
        }
    }
    modTwoMatrix(CH,4,8);
    int syndrome = syndromeLookup(p1->H, CH);

    bitFlip(cHat[0],8,syndrome);
    int (*m1)[4] = malloc(sizeof(int[1][4]));
    for (int i = 0; i < 4; ++i) {
        m1[0][i] = abs(cHat[0][i] % 2);
    }
    int(*inv)[4] = malloc(sizeof(int[4][4]));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j <4 ; ++j) {
            inv[i][j] = p1->S[i][j];
        }
    }
    int(*inversedS)[4] = malloc(sizeof(int[4][4]));
    inverse(inv, inversedS);
    printf("S before \n");

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            printf("%d ,", inv[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            inversedS[i][j] = abs(inversedS[i][j]);
        }
    }
    printf("inv \n");

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            printf("%d ,", inversedS[i][j]);
        }
        printf("\n");
    }

    n_G =1;
    p_G= 4;
    multiplyMatrix(res, m1, inversedS, 1, 4, 4);
    free(CH);
    free(inv);
    free(m1);
    free(cHat);
    free(inversedS);
    free(inversedP);
    free(transpose);
}

void decrpytFull(privateKEYPTR p1)
{
    int countFOR = 0;
    int countFOR2 = 0;

    int* m = (int*) malloc(sizeof(int*)*length*2);
    FILE *fp = fopen("Chipher.txt", "r");
    if (fp == NULL)
    {
        printf("Error opening the file %s", "Chipher.txt");
        return;
    }
    // write to the text file
    printf("m \n");
    for (int i = 0; i<length*2; i++)
    {
       m[i] = fgetc( fp)-'0';
        printf("%d, ",m[i]);
    }
    printf("\n");
    // close the file
    fclose(fp);
    int(*SubArr)[8] = malloc(sizeof(int[1][8]));
    int* decodedI = (int*) malloc(sizeof(int)*length);
    for (int i = 0; i < length; ++i) {
        decodedI[i] =0;
    }

    for (int i = 0; i < length / 4; i++)
    {
        printf("subarr \n");
        for (int j = 0;j<8;j++)
        {
            SubArr[0][j] = m[countFOR2];
            countFOR2++;
            printf("%d,", SubArr[0][j]);
        }
        printf("\n");
        int(*dec)[4] = malloc(sizeof(int[1][4]));
        decrypt(p1,SubArr,dec);
        for (int j = 0; j < 4; ++j) {
            decodedI[countFOR] = abs((dec[0][j]) % 2);
            countFOR++;

        }
        printf("\n");
        free(dec);
    }

    printf("\n");
    for (int i = 0; i < length; ++i) {
        printf("%d, ",decodedI[i]);
    }
    free(m);
    free(decodedI);
    free(SubArr);
}

void encrypt(int(*c)[8],int(*m)[],struct PrivateKey* p1) {
    // Generate a random error vector
    srand(time(NULL));
    int z = rand() % 7 + 1;
    // Multiply m by GPrime

    gprime(p1);
    n_G = 1;
    p_G = 8;
    multiplyMatrix(c,m,p1->GPRIME,1,4,8);
    // Apply bit flip to the result
    modTwoMatrix(c,1,8);
    bitFlip((int *)c, 8 , z);
    // Print the encrypted matrix (c)

}
void encryptFull(int cipheri[],privateKey1* p1,int* m)
{
    int (*SubArr)[4] = malloc(sizeof(int[1][4]));;

    for(int i = 0;i<length/4;i++)
    {
        for (int j = 0;j<4;j++)
        {
            SubArr[0][j] = m[i*4+j];
        }
        int (*c)[8] = malloc(sizeof(int[1][8]));
        encrypt(c,SubArr,p1);

        for (int j = 0; j < 8; ++j) {
            cipheri[i*8+j] = c[0][j];
        }
        free(c);
    }
    free(SubArr);
    FILE *fp = fopen("Chipher.txt", "w");
    if (fp == NULL)
    {
        printf("Error opening the file %s", "Chipher.txt");
        return;
    }
    // write to the text file
    for (int i = 0; i<length*2; i++)
    {
        fputc(cipheri[i]+'0', fp);

    }

    // close the file
    fclose(fp);
}
int main() {
    char filePath[MAX_PATH]; // Buffer to store the selected file path
    int* arrayFile;
    if (chooseFile(filePath)) {
        printf("File chosen: %s\n", filePath);
        arrayFile = readFileAndConvertToBits(filePath);
        for (int i = 0; i < length; ++i) {
            printf("%d, ",arrayFile[i]);
        }
        printf("\n");

    }
    else {
        printf("User canceled the file selection.\n");
    }
    privateKEYPTR p1 = (struct PrivateKey*)malloc(sizeof(struct PrivateKey));
    int G[4][8] = {
            {1, 0, 0, 0, 0, 1, 1, 1},
            {0, 1, 0, 0, 1, 0, 1, 1},
            {0, 0, 1, 0, 1, 1, 0, 1},
            {0, 0, 0, 1, 1, 1, 1, 0}
    };

    int H[ROWS_H][COLS_H] = {
            {0, 1, 1, 1, 1, 0, 0, 0},
            {1, 0, 1, 1, 0, 1, 0, 0},
            {1, 1, 0, 1, 0, 0, 1, 0},
            {1, 1, 1, 0, 0, 0, 0, 1}
    };
    p1->P  = malloc(sizeof(int[8][8]));
    p1->G  = malloc(sizeof(int[4][8]));
    p1->H  = malloc(sizeof(int[4][8]));
    p1->S  = malloc(sizeof(int[4][4]));
    p1->GPRIME  = malloc(sizeof(int[4][8]));
    for (int t = 0; t < ROWS_H; t++) {
        for (int j = 0; j < COLS_H; j++) {
            (p1->H)[t][j] = H[t][j];
        }
    }
    for (int t = 0; t < ROWS_H; t++) {
        for (int j = 0; j < COLS_H; j++) {
            (p1->G)[t][j] = G[t][j];
        }
    }
    /*
    int G[1][4] = {{1, 0, 0, 0}};
    int G_1[1][4] = {{0, 1, 1, 1}};
    int G_2[1][4] = {{0, 1, 0, 0}};
    int G_3[1][4] = {{ 1, 0, 1, 1}};
    int G_4[1][4] = {{0, 0, 1, 0}};
    int G_5[1][4] = {{1, 1, 0, 1}};
    int G_6[1][4] = {{0, 0, 0, 1}};
  //  int G_7[1][4] = {{ 1, 1, 1, 0}};
    // Assign the matrices to the private key
    int i=0;
    for (i = 0; i < ROWS_G; i++) {
        for (int j = 0; j < COLS_G; j++) {
            if(i==0&&j<4) {privateKey.G[i][j] = G[0][j];}
            if(i==0&&j>3) {privateKey.G[i][j] = G_1[0][j-4];}
            if(i==1&&j<4) {privateKey.G[i][j] = G_2[0][j];}
            if(i==1&&j>3) {privateKey.G[i][j] = G_3[0][j-4];}
            if(i==2&&j<4) {privateKey.G[i][j] = G_4[0][j];}
            if(i==2&&j>3) {privateKey.G[i][j] = G_5[0][j-4];}
            if(i==3&&j<4) {privateKey.G[i][j] = G_6[0][j];}
            if(i==3&&j>3&&j<6) {privateKey.G[i][j] = G_5[0][j-4];}
            if(i==3&&j==6) {privateKey.G[i][j] = 1;}
            if(i==3&&j==7) {privateKey.G[i][j] = 0;}


        }
    }
     */
    /* for (i = 0; i < 4; i++) {
         for (int j = 0; j < 8; j++) {
             printf("%d, ", privateKey.G[i][j]);
         }
         printf("\n");
     } */


    printf("h worked");
    genSMatrix(p1,4);
    genPMatrix(p1,8,0);
    printf("\n show p matrix \n");
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            printf("%d,", p1->P[i][j]);
        }
        printf("\n");
    }
    printf("llen is : %d \n",length);
    int* cipheri = (int*)malloc((sizeof(int)*length*2));
    encryptFull(cipheri,p1,arrayFile);
    for (int i = 0; i < length*2; ++i) {
        printf("%d, ",cipheri[i]);
    }

    printf("\n");
    decrpytFull(p1);
    free(cipheri);
    free(arrayFile);
    free(p1->G);
    free(p1->H);
    free(p1->S);
    free(p1->P);
    free(p1->GPRIME);
    free(p1);
    return 0;
}

