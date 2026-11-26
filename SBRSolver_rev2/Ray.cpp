#include "Ray.h"

Paths::Paths()
{

}

Paths::~Paths()
{
	for (auto i = rayPaths.begin(); i != rayPaths.end(); ++i)
	{
		i->clear();
	}
	rayPaths.clear();
}
