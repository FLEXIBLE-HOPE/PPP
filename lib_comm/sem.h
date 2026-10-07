/*
 * sem.h
 *
 *  Created on: Jul 16, 2013
 *      Author: jguo
 */
#ifndef SEM_H_
#define SEM_H_


void sem_free(int semid, int semnum);
void sem_lock(int semid, int semnum);
int sem_init_l(int nsems);

#endif /* SEM_H_ */
