#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gmp.h>

#define A_COEFFICENT 2147483647

#define MOD 2147483641
#define WIDTH 3
#define Curve448_Param 156326
struct TreeNode
{
    bool passed;
    struct TreeNode *left;
    struct TreeNode *right;
};
typedef struct
{
    mpz_t x;
    mpz_t z;

} Point;

mpz_t zero;
mpz_t one;
mpz_t two;
mpz_t three;
mpz_t four;
mpz_t mod;
mpz_t modMinusTwo;
mpz_t curve_param;
mpz_t curve_param_plus_two;

mpz_t TwoXcurve_param;
mpz_t A_edited_param;
// Structure to represent a point on the elliptic curve

Point zero_point;

void initializeMyStruct(Point *myStruct)
{
    mpz_init(myStruct->x);
    mpz_init(myStruct->z);
}
struct TreeNode *createNode()
{
    struct TreeNode *newNode = (struct TreeNode *)malloc(sizeof(struct TreeNode));
    if (newNode == NULL)
    {
        perror("Memory allocation failed");
        exit(1);
    }
    newNode->passed = true;
    newNode->left = NULL;
    newNode->right = NULL;
    return newNode;
}
bool compareArrays(bool arr1[], bool arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false; // Arrays are not equal
        }
    }
    return true; // Arrays are equal
}
void traverse(struct TreeNode *node, bool *array, int index)
{
    if (node->left != NULL)
    {
        // printf("%s\n", "left");
        array[index] = false;
        index = index + 1;
        node = node->left;
        traverse(node, array, index);
    }
    else if (node->right != NULL)
    {
        // printf("%s\n", "right");
        array[index] = true;
        index = index + 1;
        node = node->right;
        traverse(node, array, index);
    }
    else
    {
        // printf("End of the way\n");
        return;
    }
}

int *decimalToBinary(int decimal, int *scalar_bits)
{
    int i = 0;
    while (decimal > 0)
    {
        scalar_bits[i] = decimal % 2;
        decimal /= 2;
        i++;
    }
    return scalar_bits;
}

void modular_exponentiation(mpz_t ans, mpz_t a, mpz_t b)
{
    mpz_set_str(ans, "1", 10);
    mpz_t temp;
    mpz_init(temp);

    mpz_mod(a, a, mod);
    while (mpz_cmp(b, zero) > 0)
    {
        mpz_mod(temp, b, two);
        if (mpz_cmp(temp, one) == 0)
        { // b%2 == 1
            mpz_mul(temp, ans, a);
            mpz_mod(ans, temp, mod); //  (ans * a) % MOD;
        }
        mpz_mul(temp, a, a);
        mpz_mod(a, temp, mod); // a = (a * a) % MOD;
        mpz_fdiv_q(b, b, two);
    }
    mpz_clear(temp);
}

void inverse(mpz_t ans, mpz_t a)
{
    modular_exponentiation(ans, a, modMinusTwo);
}
int CC(struct TreeNode *root_1, struct TreeNode *root_2)
{
    bool traverese_1[32];
    bool traverese_2[32];
    traverse(root_1, traverese_1, 0);
    traverse(root_2, traverese_2, 0);
    return compareArrays(traverese_1,traverese_2,2);
}

struct TreeNode *error_detection(int scalar_bits_length, int *scalar_bits)
{

    for (int i = 0; i < scalar_bits_length; i++)
    {
        struct TreeNode *root = createNode();
        struct TreeNode *now = root;
        if (scalar_bits[scalar_bits_length] == 0)
        {
            now->left = createNode();
            now = now->left;
        }
        else
        {
            now->right = createNode();
            now = now->right;
        }

        return root;
    }
}

Point point_add_Montgomery(Point m, Point n, Point p)
{
    if (!mpz_cmp(m.x, zero_point.x) && !mpz_cmp(m.z, zero_point.z))
    {
        return n;
    }
    if (!mpz_cmp(n.x, zero_point.x) && !mpz_cmp(n.z, zero_point.z))
    {
        return m;
    }
    Point result;
    initializeMyStruct(&result);
    mpz_t temp;
    mpz_init(temp);
    mpz_t temp2;
    mpz_init(temp2);
    mpz_t temp3;
    mpz_init(temp3);

    mpz_sub(temp, m.x, m.z);    // x_m - z_m
    mpz_add(temp2, n.x, n.z);   // x_n + z_n
    mpz_mul(temp, temp, temp2); // (x_m - z_m) * (x_n + z_n)

    mpz_add(temp2, m.x, m.z);     // x_m + z_m
    mpz_sub(temp3, n.x, n.z);     // x_n - z_n
    mpz_mul(temp2, temp2, temp3); // ( x_m + z_m) * (x_n - z_n)

    mpz_add(temp, temp, temp2);  // (x_m - z_m) * (x_n + z_n) + (x_m + z_m) * (x_n - z_n)
    mpz_add(temp2, temp, temp2); // (x_m - z_m) * (x_n + z_n) - (x_m + z_m) * (x_n - z_n)

    mpz_mul(temp, temp, temp);
    mpz_mul(temp2, temp2, temp2);

    mpz_mul(temp, temp, p.z);
    mpz_mul(temp2, temp2, p.x);

    mpz_mod(result.x, temp, mod);
    mpz_mod(result.z, temp2, mod);
    return result;
}

Point point_doubling_Montgomery(Point n)
{
    Point result;
    initializeMyStruct(&result);
    mpz_t temp;
    mpz_init(temp);
    mpz_t temp2;
    mpz_init(temp2);
    mpz_t temp3;
    mpz_init(temp3);
    mpz_t temp4;
    mpz_t A_2;
    mpz_init(A_2);
    mpz_set_str(A_2, "156328", 10); // A+2
    mpz_init(temp4);
    mpz_add(temp, n.x, n.z);      // x_n + z_n
    mpz_mul(temp, temp, temp);    // (x_n + z_n)^2
    mpz_sub(temp2, n.x, n.z);     // (x_n - z_n)
    mpz_mul(temp2, temp2, temp2); // (x_n - z_n)^2
    mpz_sub(temp3, temp, temp2);  // 4*x_n*z_n

    mpz_mul(temp4, temp, temp2); // (x_n + z_n)^2 * (x_n - z_n)^2
    mpz_mod(result.x, temp4, mod);

    mpz_mul(temp, A_2, temp3);  //(A+2) * 4*x_n*z_n
    mpz_div(temp, temp, four);  //(A+2) * 4*x_n*z_n / 4
    mpz_add(temp, temp2, temp); //((A+2) * 4*x_n*z_n / 4) + (x_n - z_n)^2
    mpz_mul(temp, temp, temp3); //(((A+2) * 4*x_n*z_n / 4) + (x_n - z_n)^2)*(4*x_n*z_n)
    mpz_mod(result.z, temp, mod);

    return result;
}

Point Montgomery_ladder(Point point, int scalar)
{

    Point R0;
    initializeMyStruct(&R0);
    R0 = zero_point;
    Point R1;
    initializeMyStruct(&R1);
    R1 = point;
    int scalar_bits[32];
    int scalar_bits_length = log2(scalar) + 1;
    decimalToBinary(scalar, scalar_bits);
    for (int i = 1; i <= scalar_bits_length; i++)
    {
        if (scalar_bits[scalar_bits_length - i] == 0)
        {
            R1 = point_add_Montgomery(R0, R1, point);
            R0 = point_doubling_Montgomery(R0);
        }
        else
        {
            R0 = point_add_Montgomery(R0, R1, point);
            R1 = point_doubling_Montgomery(R1);
        }
    }

    return R0;
}

Point Montgomery_ladder_Error_detection(Point point, int scalar)
{

    Point R0;
    initializeMyStruct(&R0);
    R0 = zero_point;
    Point R1;
    initializeMyStruct(&R1);
    R1 = point;
    int scalar_bits[32];
    int scalar_bits_length = log2(scalar) + 1;
    decimalToBinary(scalar, scalar_bits);
    struct TreeNode *before_operation_root = error_detection(scalar_bits_length,scalar_bits);
    struct TreeNode *during_operation_root = createNode();
    struct TreeNode *now = during_operation_root;
    for (int i = 1; i <= scalar_bits_length; i++)
    {
        if (scalar_bits[scalar_bits_length - i] == 0)
        {
            now->left = createNode();
            now = now->left;
            R1 = point_add_Montgomery(R0, R1, point);
            R0 = point_doubling_Montgomery(R0);
        }
        else
        {
            now->right = createNode();
            now = now->right;
            R0 = point_add_Montgomery(R0, R1, point);
            R1 = point_doubling_Montgomery(R1);
        }
    }
    bool error = CC(before_operation_root, during_operation_root);
    // if(error == false){
    //     printf("Error happened\n");
    // }else{
    //     printf("No Error happened\n");
    // }
    return R0;
}

int main()
{
    mpz_init(zero);
    mpz_set_str(zero, "0", 10);
    mpz_init(one);
    mpz_set_str(one, "1", 10);
    mpz_init(two);
    mpz_set_str(two, "2", 10);
    mpz_init(three);
    mpz_set_str(three, "3", 10);
    mpz_init(four);
    mpz_set_str(four, "4", 10);
    mpz_init(mod);
    mpz_set_str(mod, "726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018365439", 10);
    initializeMyStruct(&zero_point);
    mpz_init_set_str(zero_point.x, "0", 10);
    mpz_init_set_str(zero_point.z, "0", 10);

    const unsigned char xx[1024] = "129988354940085515606626447021989972402497793038134173376702831301860351400271948975937626416854636844812623764974423336381552530833486 ";
    const unsigned char zz[1024] = "1";
    Point base_point;
    initializeMyStruct(&base_point);
    mpz_init_set_str(base_point.x, xx, 10);
    mpz_init_set_str(base_point.z, zz, 10);
    Point result;
    initializeMyStruct(&result);
    // result = point_doubling_Montgomery(base_point);
    // window_method(base_point,10);
    result = Montgomery_ladder(base_point, 15);
    gmp_printf("x is %Zd\n", result.x);
    gmp_printf("z is %Zd\n", result.z);
    result = Montgomery_ladder_Error_detection(base_point, 12);
    gmp_printf("x is %Zd\n", result.x);
    gmp_printf("z is %Zd\n", result.z);
    // point_add_Montgomery(base_point,base_point);
    // int scalar = 3; // Example scalar value
    // // Point result1 = window_method(base_point, scalar);
    // const unsigned char xx[1024] = "447468670702497611091631718201233577337781298833737865519591906877772691318306383650714189880688777691849432859684522905963692190064648";
    // const unsigned char yy[1024] = "474096517484922702670876688179267642901082320417012481397591077081035945699214874339205418130419484059898983650129747474170729891708658";
    // mpz_t x;
    // mpz_init_set_str(x,xx,10);
    // mpz_t y;
    // mpz_init_set_str(y,yy,10);

    // Point base_point = {x,y}; // Example base point on the elliptic curve
    // point_add_Montgomery(base_point,base_point);

    // mpz_t ans;
    // mpz_init(ans);
    // inverse(ans,f);
    // gmp_printf("answ is %Zd\n", ans);

    // uint128_t t = 447468670702497611091631718201233577337781298833737865519591906877772691318306383650714189880688777691849432859684522905963692190064648;
    // printf("%lld",t);
    // printf("result: (%d, %d)\n", result1.x, result1.y);

    // Point result2 = error_detection_window_method(base_point, scalar);
    // printf("error detection result: (%d, %d)\n", result2.x, result2.y);
    int scalar = 100001;

    PAPI_hl_region_begin("Montgomery_ladder");
    for(int i = 0;i<100000;i++){
        Montgomery_ladder(base_point, scalar+i);
    }
    PAPI_hl_region_end("Montgomery_ladder");

    PAPI_hl_region_begin("Montgomery_ladder_Error_detection_Tree");
    for(int i = 0;i<100000;i++){
        Montgomery_ladder_Error_detection(base_point, scalar+i);
    }
    PAPI_hl_region_end("Montgomery_ladder_Error_detection_Tree");

    mpz_clear(zero);
    mpz_clear(one);
    mpz_clear(two);
    mpz_clear(three);
    mpz_clear(mod);
    mpz_clear(modMinusTwo);

    return 0;
}
