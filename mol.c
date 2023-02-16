#include "mol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Setting up the atom with the values provided in the arguments
void atomset(atom* atom, char element[3], double* x, double* y, double* z) {
	atom->x = *x;
	atom->y = *y;
	atom->z = *z;
	int len = strlen(element);
	for (int i = 0; i < len; i++) {
		atom->element[i] = element[i];
	}
	atom->element[len] = '\0';
}

// Getting the values from the atom 
// and storing them in the pointers provided in the arguments 
void atomget(atom* atom, char element[3], double* x, double* y, double* z) {
	*x = atom->x;
	*y = atom->y;
	*z = atom->z;
	int len = strlen(atom->element);
	for (int i = 0; i < len; i++) {
		element[i] = atom->element[i];
	}
	element[len] = '\0';
}

// Setting up the bond with the values provided in the arguments
void bondset(bond* bond, atom* a1, atom* a2, unsigned char epairs) {
	bond->a1 = a1;
	bond->a2 = a2;
	bond->epairs = epairs;
}

// Getting the values from the bond 
// and storing them in the pointers provided in the arguments
void bondget(bond* bond, atom** a1, atom** a2, unsigned char* epairs) {
	*a1 = bond->a1;
	*a2 = bond->a2;
	*epairs = bond->epairs;
}

// Assigns initial memory for the molecule 
// according to the max values provided 
// and returns a pointer to it
molecule* molmalloc(unsigned short atom_max, unsigned short bond_max) {
	molecule* temp = malloc(sizeof(molecule));
	temp->atom_max = atom_max;
	temp->bond_max = bond_max;
	temp->atom_no = temp->bond_no = 0;

	temp->atoms = malloc(sizeof(atom) * atom_max);
	temp->atom_ptrs = malloc(sizeof(atom*) * atom_max);
	for (int i = 0; i < atom_max; i++) {
		temp->atom_ptrs[i] = malloc(sizeof(atom));
	}

	temp->bonds = malloc(sizeof(bond) * bond_max);
	temp->bond_ptrs = malloc(sizeof(bond*) * bond_max);
	for (int i = 0; i < bond_max; i++) {
		temp->bond_ptrs[i] = malloc(sizeof(bond));
	}
	return temp;
}

// Copies the molecule provided in the argument 
// and returns a pointer to the copy
molecule* molcopy(molecule* src) {
	molecule* dest = molmalloc(src->atom_max, src->bond_max);
	for (int i = 0; i < src->atom_no; i++) {
		molappend_atom(dest, src->atom_ptrs[i]);
	}
	for (int i = 0; i < src->bond_no; i++) {
		molappend_bond(dest, src->bond_ptrs[i]);
	}

	return dest;
}

// Frees the memory allocated for the molecule
void molfree(molecule* ptr) {
	free(ptr->atoms);
	free(ptr->atom_ptrs);
	free(ptr->bonds);
	free(ptr->bond_ptrs);
}

// Appends the atom provided in the argument 
// to the molecule provided in the argument
// and increases the atom_no of the molecule
void molappend_atom(molecule* molecule, atom* a) {
	if (molecule->atom_max == molecule->atom_no) {	// Checking if the atom_max is equal to the atom_no
		if (molecule->atom_max == 0) {
			molecule->atom_max++;		// incrementing the atom_max by 1 if it is 0
		}
		else {
			molecule->atom_max *= 2;	// Double the atom_max if it is not 0
		}

		// Reallocating the memory for the atoms and atom_ptrs
		molecule->atoms = realloc(molecule->atoms, molecule->atom_max * sizeof(atom));
		molecule->atom_ptrs = realloc(molecule->atom_ptrs, molecule->atom_max * sizeof(atom*));
		for (int i = 0; i < molecule->atom_max; i++) {
			molecule->atom_ptrs[i] = malloc(sizeof(atom));
		}
		for (int i = 0; i < molecule->atom_no; i++) {
			molecule->atom_ptrs[i] = &molecule->atoms[i];
		}
	}
	// Setting atom values and storing the reference in atom_ptrs as well
	atomset(&molecule->atoms[molecule->atom_no], a->element, &a->x, &a->y, &a->z);
	molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[molecule->atom_no];
	molecule->atom_no++;
}

// Appends the bond provided in the argument 
// to the molecule provided in the argument
// and increases the bond_no of the molecule
void molappend_bond(molecule* molecule, bond* b) {
	if (molecule->bond_max == molecule->bond_no) { // Checking if the bond_max is equal to the bond_no
		if (molecule->bond_max == 0) {
			molecule->bond_max++;	// incrementing the bond_max by 1 if it is 0
		}
		else {
			molecule->bond_max *= 2; // Double the bond_max if it is not 0
		}

		// Reallocating the memory for the bonds and bond_ptrs
		molecule->bonds = realloc(molecule->bonds, molecule->bond_max * sizeof(atom));
		molecule->bond_ptrs = realloc(molecule->bond_ptrs, sizeof(bond*) * molecule->bond_max);
		for (int i = 0; i < molecule->bond_max; i++) {
			molecule->bond_ptrs[i] = malloc(sizeof(bond));
		}
		for (int i = 0; i < molecule->bond_no; i++) {
			molecule->bond_ptrs[i] = &molecule->bonds[i];
		}
	}
	// Setting bond values and storing the reference in bond_ptrs as well
	bondset(&molecule->bonds[molecule->bond_no], b->a1, b->a2, b->epairs);
	molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[molecule->bond_no];
	molecule->bond_no++;
}


// Sorts the atoms and bonds of the molecule 
// in ascending order of z coordinate
// and returns a pointer to the sorted molecule
void molsort(molecule* molecule) {
	// Sorting the atoms in ascending order of z coordinate
	for (int i = 0; i < molecule->atom_no; i++) {
		for (int j = 0; j < molecule->atom_no - i - 1; j++) {
			if (molecule->atom_ptrs[j]->z > molecule->atom_ptrs[j + 1]->z) {
				atom* temp = molecule->atom_ptrs[j];
				molecule->atom_ptrs[j] = molecule->atom_ptrs[j + 1];
				molecule->atom_ptrs[j + 1] = temp;
			}
		}
	}

	// Sorting the bonds in ascending order of z coordinate

	// P.S: Bond's z coordinate is the average of the 
	// z coordinates of the atoms the bond is between
	for (int i = 0; i < molecule->bond_no; i++) {
		for (int j = 0; j < molecule->bond_no - i - 1; j++) {
			double avg1 = (molecule->bond_ptrs[j]->a1->z + molecule->bond_ptrs[j]->a2->z) / 2;
			double avg2 = (molecule->bond_ptrs[j + 1]->a1->z + molecule->bond_ptrs[j + 1]->a2->z) / 2;
			if (avg1 > avg2) {
				bond* temp = molecule->bond_ptrs[j];
				molecule->bond_ptrs[j] = molecule->bond_ptrs[j + 1];
				molecule->bond_ptrs[j + 1] = temp;
			}
		}
	}

}

void xrotation(xform_matrix matrix, double deg) {}

void yrotation(xform_matrix matrix, double deg) {}

void zrotation(xform_matrix matrix, double deg) {}

void mol_xform(molecule* molecule, xform_matrix matrix) {}
