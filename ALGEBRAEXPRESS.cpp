#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#define VARIABLE_RANGE 10
#define MAX_EXPRESSION_LENGTH 1000
#define MAX_NODES 10000
#define MAX_ACCUM 10000000
 
 
typedef struct Count {
	int value;
	int count;
} Count;
 
typedef struct Node Node;
 
struct Node {
	char operator;
	Node* left;
	Node* right;
	int vars;
	int max_count[VARIABLE_RANGE];
	Count** count;
};
 
 
int compare(const void* lhs, const void* rhs) {
	
	const Count* left = lhs;
	const Count* right = rhs;
	
	return left->value < right->value ? -1 :
		left->value > right->value ? 1 : 0;
}
 
Node* create_node(char operator, Node* left, Node* right,
	const int* max_count, Node** node, size_t* size) {
	
	Node* result;
	size_t i;
	size_t n = 1;
	
	for (i = 0; i < *size; ++i) {
		
		if (node[i]->operator == operator &&
			node[i]->left == left &&
			node[i]->right == right) {
			
			return node[i];
		}
	}
	
	result = malloc(sizeof(Node));
	
	if (!result) {
		exit(EXIT_FAILURE);
	}
	
	result->operator = operator;
	result->left = left;
	result->right = right;
	result->vars = (operator == 'x') ? 1 : left->vars+right->vars;
	
	for (i = 0; i < VARIABLE_RANGE; ++i) {
		result->max_count[i] =
			max_count[i] < result->vars ? max_count[i] : result->vars;
		n *= result->max_count[i]+1;
	}
	
	result->count = malloc(n*sizeof(Count*));
	
	if (!result->count) {
		exit(EXIT_FAILURE);
	}
	
	for (i = 0; i < n; ++i) {
		result->count[i] = NULL;
	}
	
	node[(*size)++] = result;
	return result;
}
 
Node* parse(const char** expression,
	const int* max_count, Node** node, size_t* size) {
	
	char last_operator[2] = {'\0'};
	Node* last_node[2] = {NULL};
	
	while (**expression && **expression != ')') {
		
		Node* current = NULL;
		
		switch (**expression) {
			
			case '+':
			case '-':
				
				if (last_node[1]) {
					last_node[0] = create_node(
						last_operator[0],last_node[0],last_node[1],
						max_count,node,size);
					last_node[1] = NULL;
				}
				
				last_operator[0] = **expression;
				break;
				
			case '*':
				
				if (last_node[1]) {
					last_operator[1] = **expression;
				} else {
					last_operator[0] = **expression;
				}
				
				break;
				
			case '(':
				
				++*expression;
				current = parse(expression,max_count,node,size);
				
			case 'x':
				
				if (!current) {
					current = create_node('x',NULL,NULL,max_count,node,size);
				}
				
				if (last_node[1]) {
					current = create_node(
						last_operator[1],last_node[1],current,
						max_count,node,size);
					last_node[1] = NULL;
				}
				
				if (last_node[0]) {
					if (last_operator[0] == '*') {
						last_node[0] = create_node(
							last_operator[0],last_node[0],current,
							max_count,node,size);
					} else {
						last_node[1] = current;
					}
				} else {
					last_node[0] = current;
				}
				
				break;
		}
		
		++*expression;
	}
	
	if (last_node[1]) {
		last_node[0] = create_node(
			last_operator[0],last_node[0],last_node[1],
			max_count,node,size);
	}
	
	return last_node[0];
}
 
void merge(char operator,
	const Count* left, const Count* right, Count** accum_ptr) {
	
	const Count* current_left = left;
	const Count* current_right;
	
	switch (operator) {
		
		case '+':
			
			while (current_left->count) {
				
				current_right = right;
				
				while (current_right->count) {
					
					(*accum_ptr)->value =
						current_left->value+current_right->value;
					(*accum_ptr)->count =
						current_left->count*current_right->count;
					++*accum_ptr;
					
					++current_right;
				}
				
				++current_left;
			}
			
			break;
		
		case '-':
			
			while (current_left->count) {
				
				current_right = right;
				
				while (current_right->count) {
					
					(*accum_ptr)->value =
						current_left->value-current_right->value;
					(*accum_ptr)->count =
						current_left->count*current_right->count;
					++*accum_ptr;
					
					++current_right;
				}
				
				++current_left;
			}
			
			break;
		
		case '*':
			
			while (current_left->count) {
				
				current_right = right;
				
				while (current_right->count) {
					
					(*accum_ptr)->value =
						current_left->value*current_right->value;
					(*accum_ptr)->count =
						current_left->count*current_right->count;
					++*accum_ptr;
					
					++current_right;
				}
				
				++current_left;
			}
			
			break;
	}
}
 
void compute(Node* node, int* count, Count** count_ptr, Count* accum);
 
void accumulate(Node* node, int index, int vars_left,
	int* count_left, int* count_right,
	int offset_left, int offset_right,
	Count*** count_left_ptr, Count*** count_right_ptr,
	Count** accum_ptr) {
	
	if (vars_left > node->left->vars) {
		return;
	}
	
	if (index < VARIABLE_RANGE) {
		
		const int size = ++count_right[index];
		int i;
		
		*count_right_ptr += offset_right*size;
		
		for (i = 0; i < size; ++i) {
			
			--count_right[index];
			*count_right_ptr -= offset_right;
			
			accumulate(node,index+1,vars_left+i,count_left,count_right,
				offset_left*(node->left->max_count[index]+1),
				offset_right*(node->right->max_count[index]+1),
				count_left_ptr,count_right_ptr,accum_ptr);
				
			++count_left[index];
			*count_left_ptr += offset_left;
		}
		
		*count_left_ptr -= offset_left*size;
		count_left[index] = 0;
		count_right[index] = size-1;
		
	} else if (vars_left == node->left->vars) {
		
		compute(node->left,count_left,*count_left_ptr,*accum_ptr);
		compute(node->right,count_right,*count_right_ptr,*accum_ptr);
		
		merge(node->operator,**count_left_ptr,**count_right_ptr,accum_ptr);
	}
}
 
void compute(Node* node, int* count, Count** count_ptr, Count* accum) {
	
	if (*count_ptr) {
		return;
	}
	
	if (node->operator == 'x') {
		
		int i;
		
		for (i = 0; i < VARIABLE_RANGE; ++i) {
			
			if (count[i]) {
				
				*count_ptr = malloc(2*sizeof(Count));
				
				if (!*count_ptr) {
					exit(EXIT_FAILURE);
				}
				
				(*count_ptr)[0].value = i;
				(*count_ptr)[0].count = 1;
				(*count_ptr)[1].count = 0;
				break;
			}
		}
		
	} else {
		
		int count_left[VARIABLE_RANGE] = {0};
		int* count_right = count;
		
		Count** count_left_ptr = node->left->count;
		Count** count_right_ptr = node->right->count;
		
		Count* accum_ptr = accum;
		Count* old_ptr = accum+1;
		Count* new_ptr = accum;
		size_t size;
		
		accumulate(node,0,0,count_left,count_right,1,1,
			&count_left_ptr,&count_right_ptr,&accum_ptr);
		
		qsort(accum,accum_ptr-accum,sizeof(Count),compare);
		
		while (old_ptr < accum_ptr) {
			
			if (old_ptr->value == new_ptr->value) {
				new_ptr->count += old_ptr->count;
				++old_ptr;
			} else {
				++new_ptr;
				*new_ptr = *old_ptr;
				++old_ptr;
			}
		}
		
		size = new_ptr-accum+1;
		*count_ptr = malloc((size+1)*sizeof(Count));
		
		if (!*count_ptr) {
			exit(EXIT_FAILURE);
		}
		
		memcpy(*count_ptr,accum,size*sizeof(Count));
		(*count_ptr)[size].count = 0;
	}
}
 
void destroy(Node* node) {
	
	size_t n = 1;
	size_t i;
 
	for (i = 0; i < VARIABLE_RANGE; ++i) {
		n *= node->max_count[i]+1;
	}
 
	for (i = 0; i < n; ++i) {
		free(node->count[i]);
	}
 
	free(node->count);
	free(node);
}
 
int main() {
	
	Count* accum = malloc(MAX_ACCUM*sizeof(Count));
	int cases;
	int i;
	
	if (!accum) {
		exit(EXIT_FAILURE);
	}
	
	scanf("%d",&cases);
	
	for (i = 0; i < cases; ++i) {
		
		int vars;
		int low;
		int high;
		int count[VARIABLE_RANGE] = {0};
		char expression[MAX_EXPRESSION_LENGTH];
		const char* expression_ptr = expression;
 
		Node* root;
		Node* node[MAX_NODES];
		size_t nodes = 0;
 
		Count* result = NULL;
		Count* current;
		int sum = 0;
		int j;
		
		scanf("%d %d %d",&vars,&low,&high);
		
		for (j = 0; j < vars; ++j) {
			int x;
			scanf("%d",&x);
			++count[x];
		}
		
		gets(expression);
		gets(expression);
 
		root = parse(&expression_ptr,count,node,&nodes);
		compute(root,count,&result,accum);
		current = result;
		
		while (current->count) {
			if (current->value >= low && current->value <= high) {
				sum += current->count;
			}
			++current;
		}
		
		printf("%d\n",sum);
 
		free(result);
 
		for (j = 0; j < (int)nodes; ++j) {
			destroy(node[j]);
		}
	}
	
	return EXIT_SUCCESS;
}  