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
 *	This software is distributed under the GNU General Public License.
 *   
 * ------------------------------------------------------------------------- */



#ifndef linked_list_inc
#define linked_list_inc
//********
class node
//********
{
	friend class list;
	friend class cluster;

        private:
                int value;
                node* next;
                node* last;
        public:
                node(void)
                : next(NULL), last(NULL)
                {}

                node(int val)
                : value(val), next(NULL), last(NULL)
                {}

                node(int val, node* nex, node* las)
                : value(val), next(nex), last(las)
                {}


                int get_value(void) {return value;}
};

//********
class list
//********
{
	friend class cluster;
	
	private:
                node* front;
                node* end;
        public:
                list(void)
                : front(NULL), end(NULL)
                {}

                list(int val);
	
                void append(int val);
                void remove(int val);
                void clear();
};
#endif
