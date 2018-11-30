/*
 * genericState.h
 *
 *  Created on: 30.11.2018
 *      Author: ecos
 */

#ifndef INC_GENERICSTATE_H_
#define INC_GENERICSTATE_H_

#include <memory>

template <class State>
using StateRef = std::unique_ptr<State>;

template <typename StateMachine, class State>
class GenericState
{
public:
    explicit GenericState(StateMachine &m, StateRef<State> &state) :
        m(m), state(state) {}

    template <class ConcreteState>
    static void init(StateMachine &m, StateRef<State> &state) {
        state = StateRef<State>(new ConcreteState(m, state));
        state->entry();
    }

protected:
    template <class ConcreteState>
    void change() {
        exit();
        init<ConcreteState>(m, state);
    }

    void reenter() {
        exit();
        entry();
    }

private:
    virtual void entry() {}
    virtual void exit() {}

protected:
    StateMachine &m;

private:
    StateRef<State> &state;
};

#endif /* INC_GENERICSTATE_H_ */
