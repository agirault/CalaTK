/*
*
*  Copyright 2011, 2012 by the CALATK development team
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
*
*/

#include "CImageManager.txx"

namespace CALATK
{

template class CImageManager< float, 1 >;
template class CImageManager< float, 2 >;
template class CImageManager< float, 3 >;
template class CImageManager< double, 1 >;
template class CImageManager< double, 2 >;
template class CImageManager< double, 3 >;

} // end namespace CALATK
