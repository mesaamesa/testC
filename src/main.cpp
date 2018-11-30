/*
 * main.cpp
 *
 *  Created on: 30.11.2018
 *      Author: ecos
 */

#include "stateMachine.h"

int main()
{
    Machine m;

    m.start();

    m.bringDown();
    m.bringDown();
    m.liftUp();
    m.liftUp();
    m.turnRight();
    m.paint(Machine::BLUE);
    m.turnLeft();
    m.paint(Machine::RED);
    m.paint(Machine::BLUE);

    return 0;
}


