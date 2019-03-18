// lnx - 1 = ln(x/e)

// a^2*x + b*x + c = 0
// a = 0.001 c = 0.001 b = -20

#include <stdio.h>
#include <math.h>

struct ans{
    double x1;
    double x2;
};

struct ans unstable_solve(double a, double b, double c){
    double delta = b*b - 4*a*c;
    struct ans answer;
    if(delta < 0){
        answer.x1 = HUGE_VAL;
        answer.x2 = HUGE_VAL;
        return answer;    // inf on error
    }
    delta = sqrt(delta);
    answer.x1 = (-b+delta)/(2*a);
    answer.x2 = (-b-delta)/(2*a);
    return answer;    // same values if only 1 answer
}

struct ans stable_solve(double a, double b, double c){
    double delta = b*b - 4*a*c;
    struct ans answer;
    if(delta < 0){
        answer.x1 = HUGE_VAL;
        answer.x2 = HUGE_VAL;
        return answer;    // inf on error
    }
    delta = sqrt(delta);
    if(b<0){
        answer.x1 = (-b+delta)/(2*a);
        answer.x2 = c/a/answer.x1;
    }
    else{
        answer.x2 = (-b-delta)/(2*a);
        answer.x1 = c/a/answer.x2;
    }
    return answer;    // same values if only 1 answer
}

int main(int argc, char **argv){

    double a = 0.00001;
    double b = -2000;
    double c = 0.00001;


    printf("a=%f, b=%f, c=%f\n\n",a,b,c);
    struct ans unstable_sol = unstable_solve(a,b,c);
    printf("UNSTABLE:\n%0.18f, %0.18f\n",unstable_sol.x1,unstable_sol.x2);
    struct ans stable_sol = stable_solve(a,b,c);
    printf("STABLE:\n%0.18f, %0.18f\n",stable_sol.x1,stable_sol.x2);

    return 0;
}