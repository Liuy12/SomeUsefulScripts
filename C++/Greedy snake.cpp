#include<iostream>
#include<windows.h>
#include<conio.h>
#include<time.h>
#include<stdlib.h>
using namespace std;
#define COUT cout<<"*"
struct Body
{
	int x, y;
	Body *next;
};
int Food = 0, Food_x, Food_y, Direction = 4, Grade, Score = 0;
Body *body001 = new Body[sizeof(Body)];
void Position(int x, int y)
{
	COORD pos = { y - 1, x - 1 };
	HANDLE Out = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleCursorPosition(Out, pos);
}
void CreatSneak()
{
	Body *body002 = new Body[sizeof(Body)];
	Body *body003 = new Body[sizeof(Body)];
	Body *body004 = new Body[sizeof(Body)];
	body001->x = 11; body001->y = 35;
	body002->x = 11; body002->y = 33;
	body003->x = 11; body003->y = 31;
	body004->x = 11; body004->y = 29;
	body001->next = body002;
	body002->next = body003;
	body003->next = body004;
	body004->next = NULL;
}
void Prt_Map()
{
	int i;
	for (i = 1; i <= 80; i += 2)
	{
		Position(1, i); COUT;
		Position(24, i); COUT;
	}//x:2~23 共22个
	for (i = 1; i <= 24; i++)
	{
		Position(i, 1); COUT;
		Position(i, 57); COUT;//y:3~55共27个
		Position(i, 79); COUT;
	}
}
void Prt_Sneak()
{
	Body *p; p = body001;
	while (p != NULL)
	{
		Position(p->x, p->y); COUT;
		p = p->next;
	}
	Position(Food_x, Food_y); COUT;
	Position(7, 63); cout << "Score: " << Score;
	Position(10, 63); cout << "Grade: " << Grade;
}
int JudgeOver()
{
	Body *p;
	p = body001->next;
	if (body001->y == 1 || body001->y == 57 || body001->x == 1 || body001->x == 24)
		return 1;
	while (!(p->x == body001->x&&p->y == body001->y))
	{
		if (p->next == NULL)return 0; p = p->next;
	}
	return 1;
}
void Creat_Food()
{
	if (Food == 0)
	{
		srand((int)time(0));
		Food_x = rand() % 21 + 2;
		int temp = rand() % 52 + 3;
		if (temp % 2 == 0)
			Food_y = temp + 1;
		else Food_y = temp;
		Food = 1;
	}
}
void Move()
{
	Body *p = new Body[sizeof(Body)], *q, *temp = new Body[sizeof(Body)];
	if (Direction == 1)
	{
		if (body001->x == (Food_x + 1) && body001->y == Food_y)
		{
			temp->x = Food_x; temp->y = Food_y; temp->next = body001; body001 = temp; Food = 0; Score += 5;
		}
		else
		{
			temp->x = body001->x - 1; temp->y = body001->y; temp->next = body001; body001 = temp;
			q = body001;
			while ((q->next)->next != NULL)q = q->next;
			Position((q->next)->x, (q->next)->y); cout << " ";
			delete(q->next); q->next = NULL;
		}
	}
	if (Direction == 2)
	{
		if (body001->x == Food_x&&body001->y == (Food_y + 2))
		{
			temp->x = Food_x; temp->y = Food_y; temp->next = body001; body001 = temp; Food = 0; Score += 5;
		}
		else
		{
			temp->x = body001->x; temp->y = body001->y - 2; temp->next = body001; body001 = temp;
			q = body001;
			while ((q->next)->next != NULL)q = q->next;
			Position((q->next)->x, (q->next)->y); cout << " ";
			delete(q->next); q->next = NULL;
		}
	}
	if (Direction == 3)
	{
		if (body001->x == (Food_x - 1) && body001->y == Food_y)
		{
			temp->x = Food_x; temp->y = Food_y; temp->next = body001; body001 = temp; Food = 0; Score += 5;
		}
		else
		{
			temp->x = body001->x + 1; temp->y = body001->y; temp->next = body001; body001 = temp;
			q = body001;
			while ((q->next)->next != NULL)q = q->next;
			Position((q->next)->x, (q->next)->y); cout << " ";
			delete(q->next); q->next = NULL;
		}
	}
	if (Direction == 4)
	{
		if (body001->x == Food_x&&body001->y == (Food_y - 2))
		{
			temp->x = Food_x; temp->y = Food_y; temp->next = body001; body001 = temp; Food = 0; Score += 5;
		}
		else
		{
			temp->x = body001->x; temp->y = body001->y + 2; temp->next = body001; body001 = temp;
			q = body001;
			while ((q->next)->next != NULL)q = q->next;
			Position((q->next)->x, (q->next)->y); cout << " ";
			delete(q->next); q->next = NULL;
		}
	}
}
void Game()
{
	while (1)
	{
		if (JudgeOver() == 1)return;
		if (GetAsyncKeyState(VK_UP) && Direction != 3)Direction = 1;
		if (GetAsyncKeyState(VK_LEFT) && Direction != 4)Direction = 2;
		if (GetAsyncKeyState(VK_DOWN) && Direction != 1)Direction = 3;
		if (GetAsyncKeyState(VK_RIGHT) && Direction != 2)Direction = 4;
		Creat_Food();
		Move();
		Prt_Sneak();
		Sleep(550 - Grade * 50);
	}
}
void main()
{
	Position(12, 24); cout << "Plese Select Grade:[1~10]";
	cin >> Grade; if (Grade<1 || Grade>10){ cout << "Worry!" << endl; return; }
	system("cls");
	CreatSneak();
	Prt_Map();
	Prt_Sneak();
	Game();
	system("cls");
	Position(12, 35); cout << "You Lost !" << endl;
	Position(13, 31); cout << "You Got " << Score << " Scores" << endl;
	Position(24, 29); _getch();
}