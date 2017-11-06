/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *   
 *   	https://bottogroup.wordpress.com/
 *     
 *     	Chuan Gu, c.gu@qmul.ac.uk
 *     
 *     	Copyright (2016) Botto Research Group
 *          
 *  	This software is distributed under the GNU General Public License.
 * 
 * ------------------------------------------------------------------------- */




#include"head.hpp"
#include"linked_list.hpp"
//*****************
list::list(int val)
//*****************
{
        front = new node(val);
        end = front;
}


//************************
void list::append(int val)
//************************
{

        if (front == NULL)
        {
                front = end = new node(val);
        }
        else
        {
                end->next = new node(val,NULL,end);
                end = end->next;
        }
}

//************************
void list::remove(int val)
//************************
{
        node* node_temp = front;

        while(node_temp!= NULL)
        {
                if (node_temp->value == val)
                {
			if (node_temp == front and node_temp == end)
			{
				front = NULL;
                                end   = NULL;
                                delete node_temp;
			}
                        if (node_temp == front and node_temp != end)
                        {
                                front->next->last = NULL;
                                front = front->next;
                                delete node_temp;
                        }
                        else if (node_temp == end and node_temp != front)
                        {
                                end->last->next = NULL;
                                end = end->last;
                                delete node_temp;
                        }
                        else
                        {
                                node_temp->last->next = node_temp->next;
                                node_temp->next->last = node_temp->last;
                                delete node_temp;
                        }
                }

                node_temp = node_temp->next;
        }
}

//****************
void list::clear()
//****************
{
        node* node_temp = front;

        while (node_temp != NULL)
        {
                front = front->next;
                delete node_temp;
                node_temp = front;
        }

        front = NULL;
        end   = NULL;
}

